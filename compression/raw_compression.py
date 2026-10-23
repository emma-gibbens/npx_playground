"""
Compress a SpikeGLX .bin file (+ .meta) using the WavPack and zarr 

Compression method was selected from the following paper:
Buccino AP, Winter O, Bryant D, Feng D, Svoboda K, Siegle JH. 
Compression strategies for large-scale electrophysiology data. 
J Neural Eng. 2023 Sep 18;20(5). doi: 10.1088/1741-2552/acf5a4. PMID: 37651998.

Designed to be lossless, so that the original .bin file can be restored. The .meta 
is copied to the archive, but not compressed. The .bin is read in chunks, compressed, 
and written to a zarr array. The .compression_log.json records the compression parameters, 
the SHA-256 of the original .bin (if HASH_SOURCE = True), 
the environment (pip freeze, conda env, conda list --explicit) and a copy of this script, 
so that the archive can be decompressed in the future even if the environment has changed.
A validation check is performed after compression, reading the zarr archive and comparing it 
to the original .bin in 1 s chunks (first, last, and 50 random chunks).
Run in the same environment (or compatible package versions) used for decompression; 
see the 'environment' block in the log and the .pip_freeze.txt next to the archive.

Inputs: 
- bin_file: Path to the original SpikeGLX .bin file
- meta_file: Path to the corresponding .meta file
- out_zarr: Path to the output zarr archive
- compressor_type: WavPack compressor instance
- N_JOBS: Number of parallel jobs for compression
- CHUNK_DURATION: Duration of each chunk to read from the .bin file
- HASH_SOURCE: Whether to compute and store the SHA-256 hash of the original .bin file

Outputs:
- out_zarr (.zarr):The compressed zarr archive containing the traces, will be a folder
- .meta: A copy of the original .meta file for reference
- .compression_log.json: A log file containing compression parameters, environment info, and SHA-256 hash of the original .bin file (if HASH_SOURCE = True)
- .pip_freeze.txt: A text file containing the output of pip freeze for the current environment
- .conda_env.yml: A YAML file containing the conda environment (if available)
- .conda_explicit.txt: A text file containing the output of conda list --explicit for the current environment (if available)
- .compress_script.py: A copy of this script for future reference

EFBG Oct. 2026
"""

import hashlib
import json
import platform
import shutil
import socket
import subprocess
import sys
from datetime import datetime
from importlib import metadata
from pathlib import Path
import os

import numpy as np
import spikeinterface.full as si
from wavpack_numcodecs import WavPack

########################################################################################################################
# only necessary user input
bin_file = Path(r"D:\Temp\catgt_M671_2026_08_26_g0\M671_2026_08_26_g0_tcat.imec0.ap.bin")
meta_file = bin_file.with_suffix(".meta")  # the .meta sitting next to the .bin
out_zarr = Path(r"D:\Temp\zarr\catgt_M671_2026_08_26_g0_tcat.imec0.ap.zarr")  # explicit .zarr, new location
compressor_type = WavPack(level=2)
N_JOBS = 14
CHUNK_DURATION = "1s"
HASH_SOURCE = True  # sha256 of the original .bin (a few min for a big file); lets you check a future decompression byte-for-byte
########################################################################################################################

def read_meta(meta_path):
    """Parse a SpikeGLX .meta file into a dict of strings."""
    meta = {}
    with open(meta_path, "r") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if "=" in line:
                key, value = line.split("=", 1)
                meta[key.lstrip("~")] = value
    return meta


def channel_names_from_meta(meta):
    """Channel names in file order, from snsChanMap, e.g. '(AP0;0:0)(AP1;1:1)...(SY0;1536:1536)'."""
    entries = meta["snsChanMap"].split(")")[1:]  # [0] is the '(1536,0,4' header
    return [e.lstrip("(").split(";")[0] for e in entries if e.strip()]


def sha256_file(path, block=64 * 1024 * 1024):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(block):
            h.update(chunk)
    return h.hexdigest()


def pkg_version(name):
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return None


RESTORE_INSTRUCTIONS = """\
To decompress (needs wavpack_numcodecs, numcodecs, zarr versions from 'environment'
below, or a compatible release; see the .pip_freeze.txt next to this file):

    import spikeinterface.full as si
    import wavpack_numcodecs  # registers the codec with numcodecs
    rec = si.read_zarr(r"<path to .zarr>")
    n = rec.get_num_samples(); step = int(rec.sampling_frequency) * 10
    with open("restored.ap.bin", "wb") as f:
        for s in range(0, n, step):
            f.write(rec.get_traces(start_frame=s, end_frame=min(s + step, n)).tobytes())

Raw int16, shape (samples, channels), C-order, same as SpikeGLX .bin.
sha256 of restored.ap.bin should equal source.sha256 below. Compression is lossless.
"""


def write_log(log_path, log):
    with open(log_path, "w") as f:
        json.dump(log, f, indent=2, default=str)


def main():
    started = datetime.now()
    meta = read_meta(meta_file)

    n_chan = int(meta["nSavedChans"])
    fs = float(meta["imSampRate"])
    names = channel_names_from_meta(meta)
    assert len(names) == n_chan, f"{len(names)} names in snsChanMap vs nSavedChans={n_chan}"

    # AP gain (uV/bit) = Vrange / maxInt / gain; sync channels are digital, so gain 1
    n_ap, n_lf, n_sy = (int(x) for x in meta["snsApLfSy"].split(","))
    ap_gain = float(meta["imAiRangeMax"]) / float(meta["imMaxInt"]) / float(meta["imChan0apGain"]) * 1e6
    gains = np.concatenate([np.full(n_ap + n_lf, ap_gain), np.ones(n_sy)])

    # consistency checks between meta and the actual file
    n_samples, remainder = divmod(bin_file.stat().st_size, n_chan * 2)
    assert remainder == 0, "file size is not a whole number of samples for this channel count"
    print(f"{n_chan} channels, {n_samples} samples, {n_samples / fs:.1f} s, "
          f"meta says {float(meta['fileTimeSecs']):.1f} s, gain {ap_gain:.5f} uV/bit")

    # one definition of the archive path, used everywhere below
    archive = out_zarr.with_name(out_zarr.stem + ".zarr")
    log_path = archive.with_name(archive.stem + ".compression_log.json")

    log = {
        "restore_instructions": RESTORE_INSTRUCTIONS,
        "status": "compression started",
        "started": started.isoformat(timespec="seconds"),
        "machine": {"host": socket.gethostname(), "os": platform.platform(),
                    "python": sys.version},
        "source": {
            "bin_file": str(bin_file),
            "bin_size_bytes": bin_file.stat().st_size,
            "meta_file": str(meta_file),
            "sha256": sha256_file(bin_file) if HASH_SOURCE else None,
            "meta_sha256": sha256_file(meta_file),
        },
        "recording": {
            "dtype": "int16", "sampling_frequency": fs, "num_channels": n_chan,
            "num_samples": n_samples, "duration_s": n_samples / fs,
            "n_ap_lf_sy": [n_ap, n_lf, n_sy], "ap_gain_uV_per_bit": ap_gain,
            "channel_ids": names,
        },
        "compression": {
            "format": "zarr",
            "compressor_class": f"{type(compressor_type).__module__}.{type(compressor_type).__name__}",
            "compressor_config": compressor_type.get_config(),
            "save_kwargs": {"n_jobs": N_JOBS, "chunk_duration": CHUNK_DURATION},
            "output_zarr": str(archive),
        },
        "environment": {p: pkg_version(p) for p in
                        ["spikeinterface", "zarr", "numcodecs", "wavpack-numcodecs",
                         "numpy", "probeinterface", "neo"]},
    }

    # save full environment (pip) + conda environment + exact copy of this script
    archive.parent.mkdir(parents=True, exist_ok=True)
    freeze = subprocess.run([sys.executable, "-m", "pip", "freeze"],
                            capture_output=True, text=True).stdout
    archive.with_name(archive.stem + ".pip_freeze.txt").write_text(freeze)
    shutil.copy(Path(__file__), archive.with_name(archive.stem + ".compress_script.py"))
    conda_exe = os.environ.get("CONDA_EXE", "conda")
    for args, suffix in [(["env", "export"], ".conda_env.yml"),
                        (["list", "--explicit"], ".conda_explicit.txt")]:
        res = subprocess.run([conda_exe, *args, "--prefix", sys.prefix],
                            capture_output=True, text=True)
        if res.returncode == 0:
            archive.with_name(archive.stem + suffix).write_text(res.stdout)
        else:
            print(f"conda {' '.join(args)} failed: {res.stderr}")

    write_log(log_path, log)  

    rec = si.read_binary(
        bin_file,
        sampling_frequency=fs,
        dtype="int16",
        num_channels=n_chan,
        channel_ids=names,
        gain_to_uV=gains,
        offset_to_uV=np.zeros(n_chan),
    )

    rec.save(
        format="zarr",
        folder=out_zarr,
        compressor=compressor_type,
        n_jobs=N_JOBS,
        chunk_duration=CHUNK_DURATION,
    )

    # keep the meta file with the archive
    shutil.copy(meta_file, archive.with_suffix(".meta"))

    # record what zarr actually wrote (not just what I tried to do)
    out_bytes = sum(p.stat().st_size for p in archive.rglob("*") if p.is_file())
    log["compression"]["output_size_bytes"] = out_bytes
    log["compression"]["ratio"] = bin_file.stat().st_size / out_bytes
    try:
        import zarr
        arr = zarr.open(str(archive), mode="r")["traces_seg0"]
        log["compression"]["zarr_array_info"] = str(arr.info)
    except Exception as e:
        log["compression"]["zarr_array_info"] = f"could not read: {e!r}"
    log["status"] = "compressed, verification pending"
    write_log(log_path, log)

    # verify the archive against the original
    comp = si.read_zarr(archive)
    assert comp.get_num_channels() == n_chan
    assert comp.get_num_samples() == n_samples

    chunk = int(fs)  # 1 s
    rng = np.random.default_rng(0)
    starts = [0, n_samples - chunk] + list(rng.integers(0, n_samples - chunk, size=50))
    for s in starts:
        a = rec.get_traces(start_frame=int(s), end_frame=int(s) + chunk)
        b = comp.get_traces(start_frame=int(s), end_frame=int(s) + chunk)
        assert np.array_equal(a, b), f"MISMATCH in 1 s chunk starting at sample {s}"
    print(f"verification passed on {len(starts)} chunks (first, last, and random)")

    log["verification"] = {"n_chunks": len(starts), "chunk_samples": chunk,
                           "rng_seed": 0, "passed": True}
    log["status"] = "complete, verified"
    log["finished"] = datetime.now().isoformat(timespec="seconds")
    write_log(log_path, log)


if __name__ == "__main__":
    main()