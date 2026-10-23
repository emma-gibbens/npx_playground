"""
Restore a SpikeGLX .bin (+ .meta) from the WavPack/zarr archive made by raw_compression.py.

Reads the archive's .compression_log.json, decompresses the traces in chunks, hashes the
bytes as they are written, and compares against the SHA-256 of the original .bin recorded
at compression time. The file is only given its final name if size and hash check out.

Run in the same environment (or compatible package versions) used for compression;
see the 'environment' block in the log and the .pip_freeze.txt next to the archive.
"""
import hashlib
import json
import shutil
from datetime import datetime
from importlib import metadata
from pathlib import Path

import numpy as np
import spikeinterface.full as si
import wavpack_numcodecs  # noqa: F401  (importing registers the WavPack codec with numcodecs, don't delete)

########################################################################################################################
# user input 
archive = Path(r"D:\Temp\zarr\catgt_M671_2026_08_26_g0_tcat.imec0.ap.zarr")
restore_dir = Path(r"D:\Temp\uncompressed\catgt_M671_2026_08_26_g0")   # folder your pipeline reads as CATGT_OUTPUT_DIR
restored_stem = "M671_2026_08_26_g0_tcat.imec0.ap"        # ORIGINAL SpikeGLX name, without .bin/.meta
CHUNK_SECONDS = 10                                        # ~230 MB per chunk for 385 channels
OVERWRITE = False

########################################################################################################################

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


def main():
    started = datetime.now()

    # locate archive pieces written by raw_compression.py
    log_path = archive.with_name(archive.stem + ".compression_log.json")
    meta_archived = archive.with_suffix(".meta")
    for p in (archive, log_path, meta_archived):
        if not p.exists():
            raise FileNotFoundError(f"Missing archive component: {p}")

    with open(log_path, "r") as f:
        log = json.load(f)

    if log.get("status") != "complete, verified":
        print(f"WARNING: compression log status is '{log.get('status')}', not 'complete, verified'. "
              "The archive may be incomplete.")

    n_chan = log["recording"]["num_channels"]
    n_samples = log["recording"]["num_samples"]
    expected_size = log["source"]["bin_size_bytes"]
    expected_sha = log["source"].get("sha256")  # None if HASH_SOURCE was False at compression time
    expected_meta_sha = log["source"].get("meta_sha256")

    # check the archived .meta is intact 
    meta_sha = sha256_file(meta_archived)
    if expected_meta_sha and meta_sha != expected_meta_sha:
        raise RuntimeError("Archived .meta does not match the hash recorded at compression time")

    # create the output paths, check for existing files, and make the restore folder if needed
    out_bin = restore_dir / f"{restored_stem}.bin"
    out_meta = restore_dir / f"{restored_stem}.meta"
    tmp_bin = restore_dir / f"{restored_stem}.bin.partial"
    if (out_bin.exists() or out_meta.exists()) and not OVERWRITE:
        raise FileExistsError(f"{out_bin} or {out_meta} already exists (set OVERWRITE = True to replace)")
    restore_dir.mkdir(parents=True, exist_ok=True)

    # open archive and check it against the log
    rec = si.read_zarr(archive)
    assert rec.get_num_segments() == 1, "expected a single-segment recording"
    assert rec.get_num_channels() == n_chan, f"{rec.get_num_channels()} channels in zarr vs {n_chan} in log"
    assert rec.get_num_samples() == n_samples, f"{rec.get_num_samples()} samples in zarr vs {n_samples} in log"
    assert rec.get_dtype() == np.int16, f"zarr dtype is {rec.get_dtype()}, expected int16"
    assert n_samples * n_chan * 2 == expected_size, "log's sample/channel counts don't match its recorded .bin size"

    # decompress in chunks, hashing the exact bytes written
    h = hashlib.sha256()
    step = int(rec.sampling_frequency) * CHUNK_SECONDS
    n_chunks = -(-n_samples // step)
    print(f"Restoring {n_samples} samples x {n_chan} channels to {out_bin}")
    with open(tmp_bin, "wb") as f:
        for i, s in enumerate(range(0, n_samples, step), start=1):
            e = min(s + step, n_samples)
            traces = rec.get_traces(start_frame=s, end_frame=e)  # raw int16, (samples, channels)
            assert traces.dtype == np.int16
            data = np.ascontiguousarray(traces).tobytes()         # C-order, same layout as SpikeGLX .bin
            f.write(data)
            h.update(data)
            if i % 20 == 0 or i == n_chunks:
                print(f"  chunk {i}/{n_chunks}")
    restored_sha = h.hexdigest()

    # verify before giving the file its real name 
    size = tmp_bin.stat().st_size
    if size != expected_size:
        raise RuntimeError(f"Restored size {size} != original {expected_size}; left as {tmp_bin}")

    if expected_sha is None:
        hash_status = "not checked (original hash was not recorded)"
        print("WARNING: no source SHA-256 in the log, so only file size was verified.")
    elif restored_sha == expected_sha:
        hash_status = "match"
        print("SHA-256 matches the original .bin: restore is byte-identical.")
    else:
        raise RuntimeError(
            f"SHA-256 MISMATCH. restored={restored_sha} original={expected_sha}. Left as {tmp_bin}; "
            "check that package versions match the 'environment' block in the compression log."
        )

    tmp_bin.replace(out_bin)
    shutil.copy(meta_archived, out_meta)

    # record what happened, next to the archive (next to the original compression log)
    restore_log = {
        "restored_at": datetime.now().isoformat(timespec="seconds"),
        "duration_s": (datetime.now() - started).total_seconds(),
        "archive": str(archive),
        "restored_bin": str(out_bin),
        "restored_meta": str(out_meta),
        "restored_sha256": restored_sha,
        "expected_sha256": expected_sha,
        "hash_check": hash_status,
        "size_bytes": size,
        "environment": {p: pkg_version(p) for p in
                        ["spikeinterface", "zarr", "numcodecs", "wavpack-numcodecs", "numpy"]},
    }
    with open(archive.with_name(archive.stem + ".restore_log.json"), "w") as f:
        json.dump(restore_log, f, indent=2)
    print("Done.")


if __name__ == "__main__":
    main()