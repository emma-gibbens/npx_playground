"""
Saturation diagnostics for Neuropixels 2.0 (e.g. NP2020 quad-base) AP recordings.

Reads the RAW int16 data (no filtering / CMR / scaling) in chunks and produces:

  Figure 1: histogram of raw ADC codes, with zooms on both extremes.
            A true rail shows up as a pile-up of samples at the min/max code.
  Figure 2: distribution of per-chunk, per-channel max |V| (uV), the same
            distribution collapsed across channels at the quantile that matches
            `proportion`, and a time course of the fraction of channels at the rail.

It prints an estimated rail (uV) and a suggested saturation_threshold_uV.

Usage
  - Standalone: edit the parameters below and run the script.
  - From your pipeline: `from saturation_diagnostics import run_saturation_diagnostics`
    and call it with the RAW (unfiltered) SpikeInterface recording.
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# ----------------------------- parameters ------------------------------------
DATA_FOLDER = r"D:\Task4\M671\rawdata\M671_2026_08_26\M671_2026_08_26_g0\M671_2026_08_26_g0_imec0"   # folder with the .ap.bin and .ap.meta
STREAM_NAME = "imec0.ap"                # SpikeGLX stream name
CHUNK_S = 0.5            # chunk length (s) for the per-chunk max |V|
CHUNK_STRIDE = 1         # use every Nth chunk (raise to speed up; artifacts are rare, so keep 1 if you can)
START_S, END_S = 0.0, None   # restrict to a time window (s); END_S=None -> end of recording
PROPORTIONS = (0.05, 0.2, 0.5)   # proportion values to compare in Figure 2B
FALLBACK_ADC_MAX = 8192  # NP2 (14-bit) imMaxInt, used only if no .meta is found
OUT_DIR = Path("saturation_diagnostics")
# ------------------------------------------------------------------------------

CODE_OFFSET = 32768  # int16 codes -> non-negative bin index


def read_spikeglx_meta(meta_path):
    meta = {}
    with open(meta_path, "r") as f:
        for line in f:
            if "=" in line:
                k, v = line.strip().split("=", 1)
                meta[k.lstrip("~")] = v
    return meta


def find_meta(folder, stream_name):
    hits = sorted(Path(folder).rglob(f"*.{stream_name}.meta"))
    return hits[0] if hits else None


def edge_report(counts, code, inward):
    """Counts at `code` and the next few codes toward the interior (inward=+1 at the
    bottom extreme, -1 at the top), plus a ratio vs. codes 3-12 steps inside."""
    idx = code + CODE_OFFSET
    near = [int(counts[idx + inward * k]) for k in range(0, 6)]
    ref = np.mean([counts[idx + inward * k] for k in range(3, 13)])
    ratio = counts[idx] / max(ref, 1.0)
    return near, ratio


def run_saturation_diagnostics(rec, meta=None, chunk_s=CHUNK_S, chunk_stride=CHUNK_STRIDE,
                               start_s=START_S, end_s=END_S, proportions=PROPORTIONS,
                               fallback_adc_max=FALLBACK_ADC_MAX, out_dir=OUT_DIR):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fs = rec.get_sampling_frequency()
    n_total = rec.get_num_channels()
    n_ap = 1536 if n_total >= 1536 else (384 if n_total >= 384 else n_total)  # drop sync channels
    ap_ids = rec.channel_ids[:n_ap]
    gains = np.asarray(rec.get_channel_gains())[:n_ap].astype(float)   # uV per ADC code
    gain_uV = float(np.median(gains))

    if meta is not None and "imMaxInt" in meta:
        adc_max = int(meta["imMaxInt"])
        print(f"imMaxInt from meta: {adc_max}  | imAiRangeMax/Min: "
              f"{meta.get('imAiRangeMax')}/{meta.get('imAiRangeMin')}")
    else:
        adc_max = fallback_adc_max
        print(f"No meta imMaxInt available; assuming {adc_max}")
    print(f"gain_to_uV (median over AP channels): {gain_uV:.4f} uV/code "
          f"-> theoretical rail ~ +/-{adc_max * gain_uV:.0f} uV")

    n_samp = rec.get_num_samples()
    start = int(start_s * fs)
    end = n_samp if end_s is None else min(n_samp, int(end_s * fs))
    chunk = int(chunk_s * fs)
    starts = np.arange(start, end - chunk + 1, chunk * chunk_stride)
    print(f"Processing {len(starts)} chunks of {chunk_s}s over {n_ap} AP channels...")

    counts = np.zeros(65536, dtype=np.int64)
    maxabs = np.zeros((len(starts), n_ap), dtype=np.uint16)   # per chunk, per channel max |code|

    try:
        from tqdm import tqdm
        it = tqdm(enumerate(starts), total=len(starts))
    except ImportError:
        it = enumerate(starts)

    for i, s in it:
        x = rec.get_traces(start_frame=int(s), end_frame=int(s) + chunk,
                           channel_ids=ap_ids, return_scaled=False)
        xi = x.astype(np.int32)
        counts += np.bincount((xi + CODE_OFFSET).ravel(), minlength=65536)
        maxabs[i] = np.abs(xi).max(axis=0)

    # ------------------------- Figure 1: raw code histogram ---------------------
    codes = np.arange(-CODE_OFFSET, CODE_OFFSET)
    nz = codes[counts > 0]
    lo, hi = int(nz.min()), int(nz.max())
    print(f"\nObserved code range: {lo} .. {hi}   (theoretical +/-{adc_max})")

    near_hi, ratio_hi = edge_report(counts, hi, inward=-1)
    near_lo, ratio_lo = edge_report(counts, lo, inward=+1)
    print(f"Top    code {hi}: counts at hi, hi-1, ... = {near_hi}  (ratio vs interior: {ratio_hi:.1f})")
    print(f"Bottom code {lo}: counts at lo, lo+1, ... = {near_lo}  (ratio vs interior: {ratio_lo:.1f})")
    for c in (adc_max - 1, -adc_max):
        print(f"Count at theoretical extreme code {c}: {int(counts[c + CODE_OFFSET])}")

    # Heuristic only: a rail pile-up has far more samples at the extreme code than just inside it
    piled = [(c, r) for c, r, n in ((hi, ratio_hi, near_hi), (lo, ratio_lo, near_lo))
             if r > 5 and n[0] >= 100]
    if piled:
        rail_code = max(abs(c) for c, _ in piled)
        print(f"Likely pile-up at code(s) {[c for c, _ in piled]} -> rail ~ {rail_code} codes")
    else:
        rail_code = adc_max
        print("No clear pile-up at the extreme codes (this recording may not saturate, "
              "or the pile-up is at a code other than the min/max). Using theoretical rail.")
    rail_uV = rail_code * gain_uV
    suggested_uV = 0.95 * rail_uV
    print(f"\nEstimated rail: {rail_uV:.0f} uV | suggested saturation_threshold_uV (95% of rail): {suggested_uV:.0f}")

    fig1, axes = plt.subplots(1, 3, figsize=(16, 4))
    axes[0].semilogy(codes, np.maximum(counts, 0.5), lw=0.8)
    axes[0].set(title="Raw ADC code histogram", xlabel="int16 code", ylabel="samples")
    for ax, sel, title in ((axes[1], slice(hi - 60, hi + 2), f"Top extreme (max={hi})"),
                           (axes[2], slice(lo - 1, lo + 61), f"Bottom extreme (min={lo})")):
        cc = codes[(codes >= sel.start) & (codes <= sel.stop)]
        ax.semilogy(cc, np.maximum(counts[cc + CODE_OFFSET], 0.5), "o-", ms=3, lw=0.8)
        ax.set(title=title, xlabel="int16 code")
    fig1.tight_layout()
    fig1.savefig(out_dir / "fig1_code_histogram.png", dpi=150)

    # ------------------------- Figure 2: per-chunk max |V| ----------------------
    maxabs_uV = maxabs.astype(np.float64) * gains[None, :]
    bins = np.logspace(np.log10(max(gain_uV, 1.0)), np.log10(rail_uV * 1.2), 120)

    fig2, axes = plt.subplots(1, 3, figsize=(17, 4))
    axes[0].hist(maxabs_uV.ravel(), bins=bins, histtype="stepfilled", alpha=0.6)
    axes[0].set(xscale="log", yscale="log", title="Per-chunk, per-channel max |V|",
                xlabel="max |V| in chunk (uV)", ylabel="channel-chunks")

    for p in proportions:
        q = np.quantile(maxabs_uV, 1 - p, axis=1)   # value that a fraction p of channels exceed, per chunk
        axes[1].hist(q, bins=bins, histtype="step", lw=1.5, label=f"proportion={p}")
    axes[1].set(xscale="log", yscale="log", xlabel="(1 - proportion) quantile across channels (uV)",
                ylabel="chunks", title="Chunk-level view matching `proportion`")
    axes[1].legend()

    for ax in axes[:2]:
        ax.axvline(rail_uV, color="r", ls="--", label="rail")
        ax.axvline(suggested_uV, color="k", ls=":", label="suggested thr")
    axes[0].legend()

    frac_rail = (maxabs >= 0.98 * rail_code).mean(axis=1)
    axes[2].plot(starts / fs, frac_rail, lw=0.8)
    axes[2].set(title="Fraction of channels within 2% of rail, per chunk",
                xlabel="time (s)", ylabel="fraction of channels")
    fig2.tight_layout()
    fig2.savefig(out_dir / "fig2_chunk_max_distribution.png", dpi=150)

    np.savez(out_dir / "saturation_diagnostics.npz", counts=counts, maxabs_code=maxabs,
             chunk_starts=starts, gains=gains, rail_code=rail_code, fs=fs)
    print(f"\nSaved figures and arrays to {out_dir.resolve()}")
    plt.show()
    return dict(rail_uV=rail_uV, suggested_threshold_uV=suggested_uV, rail_code=rail_code,
                counts=counts, maxabs_uV=maxabs_uV)


if __name__ == "__main__":
    import spikeinterface.full as si

    rec = si.read_spikeglx(DATA_FOLDER, stream_name=STREAM_NAME)  # raw, no preprocessing
    meta_path = find_meta(DATA_FOLDER, STREAM_NAME)
    meta = read_spikeglx_meta(meta_path) if meta_path else None
    run_saturation_diagnostics(rec, meta=meta)
