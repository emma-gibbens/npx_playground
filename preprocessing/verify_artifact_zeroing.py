#!/usr/bin/env python
"""
verify_artifact_zeroing.py

Confirms that the artifact periods found by extract_lfp_and_run_ks4_update.py were really zeroed in
the data that Kilosort4 sees.

In the pipeline, the AP chain is:
    raw -> bandpass -> remove bad channels -> common reference -> silence_periods(zeros) -> save binary
so the binary saved in FAST_STORAGE_DIR\\<SESSION_ID>_imec<N>[_shank<K>]_si_preprocess is the final
data Kilosort reads. This script opens that saved recording (or, if it was cleaned up, rebuilds the same
chain without the bad-channel step) and runs three checks.

CHECK 1  Zeroing       Every sample inside every artifact block must be exactly 0 on every channel.
CHECK 2  Boundaries    Largest absolute value in a short window just before / just after each block,
                       in units of the recording's noise level. Silencing happens AFTER filtering, so
                       a detected span that is too tight (or filter ringing) would leave a big
                       transient right at the edge. This is a heuristic, reported as warnings.
CHECK 3  Full scan     (optional, slow, --full_scan) Reads the whole file and confirms that the ONLY
                       all-zero samples are the expected ones: nothing missed, nothing extra.

It also saves a figure of raw vs saved-preprocessed traces around example artifacts, so you can see
the zeroed gap with your own eyes.

Usage
-----
    python verify_artifact_zeroing.py --params D:\\path\\to\\M671_2026_08_27_lfp_ks4_params.txt --shank 0
    python verify_artifact_zeroing.py --params ...txt --shank 0 --full_scan

It reads the artifact list cached by artifacts_vs_reward.py (imec0_artifacts.csv), so run that first.
Run it once per shank you care about (--shank 0, 1, 2, 3).
"""
import argparse
import importlib.util
import logging
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

PIPELINE_FILENAME = "extract_lfp_and_run_ks4_update.py"
logger = logging.getLogger("verify_artifact_zeroing")


# ----------------------------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------------------------
def merge_blocks(starts, ends):
    """Merge overlapping/touching [start, end) sample intervals. Returns a list of [start, end]."""
    order = np.argsort(starts)
    blocks = []
    for i in order:
        s = int(starts[i])
        e = int(ends[i])
        if len(blocks) > 0 and s <= blocks[-1][1]:
            if e > blocks[-1][1]:
                blocks[-1][1] = e
        else:
            blocks.append([s, e])
    return blocks


def get_traces(rec, start, end, scaled=False, channel_ids=None):
    """Read samples [start, end). Unscaled by default so that exact zeros stay exact zeros."""
    kwargs = dict(start_frame=start, end_frame=end, return_scaled=scaled)
    if channel_ids is not None:
        kwargs["channel_ids"] = channel_ids
    if scaled:
        try:
            return rec.get_traces(**kwargs)
        except Exception:
            kwargs["return_scaled"] = False
            return rec.get_traces(**kwargs)
    return rec.get_traces(**kwargs)


def block_overlaps_window(blocks, start, end):
    for block in blocks:
        if block[0] < end and block[1] > start:
            return True
    return False


# ----------------------------------------------------------------------------------------------
# Opening the data
# ----------------------------------------------------------------------------------------------
def load_pipeline_module(path):
    spec = importlib.util.spec_from_file_location("ks4_pipeline", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def open_raw(params, probe_idx, shank_idx):
    import spikeinterface.full as si

    raw = si.read_spikeglx(params["CATGT_OUTPUT_DIR"], stream_name=f"imec{probe_idx}.ap")
    if shank_idx is not None:
        raw = raw.split_by("group")[shank_idx]
    return raw


def open_processed(pipeline, params, probe_idx, shank_idx, blocks, folder_override):
    """Open the saved preprocessed recording, or rebuild the same chain if it is gone."""
    import spikeinterface.full as si

    if shank_idx is not None:
        folder_name = f"{params['SESSION_ID']}_imec{probe_idx}_shank{shank_idx}"
    else:
        folder_name = f"{params['SESSION_ID']}_imec{probe_idx}"
    folder = folder_override or os.path.join(params["FAST_STORAGE_DIR"], folder_name + "_si_preprocess")

    if os.path.isdir(folder):
        logger.info(f"Opening saved preprocessed recording: {folder}")
        try:
            rec = si.load(folder)
        except AttributeError:
            rec = si.load_extractor(folder)
        return rec, f"saved file ({folder})"

    logger.warning(f"{folder} not found. Rebuilding the AP chain (no bad-channel removal) instead.")
    rec = open_raw(params, probe_idx, shank_idx)
    if params["AP_LOW_PASS"] is not None:
        rec = si.bandpass_filter(rec, freq_min=params["AP_HIGH_PASS"], freq_max=params["AP_LOW_PASS"])
    else:
        rec = si.highpass_filter(rec, freq_min=params["AP_HIGH_PASS"])
    if params["COMMON_REFERENCE_AP"]:
        rec, ok = pipeline.common_reference(rec, logger, params)
    periods = np.zeros(len(blocks), dtype=[("segment_index", "int64"),
                                           ("start_sample_index", "int64"),
                                           ("end_sample_index", "int64")])
    for i, block in enumerate(blocks):
        periods["start_sample_index"][i] = block[0]
        periods["end_sample_index"][i] = block[1]
    rec = si.silence_periods(rec, periods, mode="zeros")
    return rec, "REBUILT in memory (not the saved file)"


# ----------------------------------------------------------------------------------------------
# Checks
# ----------------------------------------------------------------------------------------------
def estimate_noise(rec, blocks, fs, n_chunks=8, chunk_s=1.0, seed=0):
    """Robust noise level (median(|x|)/0.6745) from random 1 s chunks that contain no artifact."""
    rng = np.random.default_rng(seed)
    n_samples = rec.get_num_samples()
    chunk = int(chunk_s * fs)
    medians = []
    attempts = 0
    while len(medians) < n_chunks and attempts < 200:
        attempts += 1
        start = int(rng.integers(0, max(1, n_samples - chunk)))
        end = start + chunk
        if block_overlaps_window(blocks, start, end):
            continue
        traces = get_traces(rec, start, end)
        medians.append(np.median(np.abs(traces)))
    if len(medians) == 0:
        raise RuntimeError("Could not find artifact-free chunks to estimate the noise level.")
    return float(np.median(medians)) / 0.6745


def check_blocks(rec, blocks, fs, pad_ms, sigma):
    """CHECK 1 (zeroing) and CHECK 2 (boundary transients) for every block."""
    n_samples = rec.get_num_samples()
    pad = int(round(pad_ms * fs / 1000.0))
    rows = []
    for number, block in enumerate(blocks):
        start = block[0]
        end = min(block[1], n_samples)
        inside = get_traces(rec, start, end)
        n_nonzero = int(np.count_nonzero(inside))
        max_inside = float(np.max(np.abs(inside))) if inside.size > 0 else 0.0

        before_start = max(0, start - pad)
        after_end = min(n_samples, end + pad)
        ratio_before = np.nan
        ratio_after = np.nan
        if start > before_start:
            before = get_traces(rec, before_start, start)
            ratio_before = float(np.max(np.abs(before))) / sigma
        if after_end > end:
            after = get_traces(rec, end, after_end)
            ratio_after = float(np.max(np.abs(after))) / sigma

        rows.append(dict(
            block=number,
            start_sample=start,
            end_sample=end,
            start_s=start / fs,
            duration_ms=(end - start) * 1000.0 / fs,
            n_values_checked=int(inside.size),
            n_nonzero_inside=n_nonzero,
            max_abs_inside=max_inside,
            edge_before_in_sigma=ratio_before,
            edge_after_in_sigma=ratio_after,
        ))
    return pd.DataFrame(rows)


def full_scan(rec, blocks, fs, chunk_s=10.0):
    """CHECK 3: the all-zero samples in the whole file must be exactly the expected ones."""
    n_samples = rec.get_num_samples()
    chunk = int(chunk_s * fs)
    expected_total = 0
    zero_total = 0
    missed_total = 0
    unexpected_total = 0
    unexpected_examples = []
    n_chunks = int(np.ceil(n_samples / chunk))
    chunk_number = 0
    for start in range(0, n_samples, chunk):
        end = min(start + chunk, n_samples)
        traces = get_traces(rec, start, end)
        all_zero = np.all(traces == 0, axis=1)

        expected = np.zeros(end - start, dtype=bool)
        for block in blocks:
            if block[0] < end and block[1] > start:
                lo = max(block[0], start) - start
                hi = min(block[1], end) - start
                expected[lo:hi] = True

        expected_total += int(expected.sum())
        zero_total += int(all_zero.sum())
        missed_total += int(np.sum(expected & ~all_zero))
        unexpected = all_zero & ~expected
        unexpected_total += int(unexpected.sum())
        if len(unexpected_examples) < 10:
            where = np.where(unexpected)[0]
            for w in where:
                if len(unexpected_examples) < 10:
                    unexpected_examples.append((start + int(w)) / fs)
        chunk_number += 1
        if chunk_number % 20 == 0:
            logger.info(f"  full scan: chunk {chunk_number}/{n_chunks}")
    return dict(expected_zero_samples=expected_total, observed_allzero_samples=zero_total,
                expected_but_not_zero=missed_total, zero_but_not_expected=unexpected_total,
                unexpected_times_s=unexpected_examples)


# ----------------------------------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------------------------------
def pick_examples(blocks, n_examples):
    """The longest block, then blocks spread evenly through the session."""
    chosen = []
    longest = 0
    for i, block in enumerate(blocks):
        if (block[1] - block[0]) > (blocks[longest][1] - blocks[longest][0]):
            longest = i
    chosen.append(longest)
    if n_examples > 1:
        positions = np.linspace(0, len(blocks) - 1, n_examples - 1)
        for p in positions:
            index = int(round(p))
            if index not in chosen:
                chosen.append(index)
    return chosen


def plot_examples(proc_rec, raw_rec, blocks, fs, example_indices, window_ms, n_channels, out_png):
    channel_ids = proc_rec.channel_ids
    picks = np.linspace(0, len(channel_ids) - 1, n_channels).astype(int)
    chosen_ids = []
    for p in picks:
        chosen_ids.append(channel_ids[p])

    half = int(round(window_ms * fs / 1000.0))
    n_rows = len(example_indices)
    fig, axes = plt.subplots(n_rows, 2, figsize=(13, 2.6 * n_rows), squeeze=False)

    for row, block_index in enumerate(example_indices):
        start, end = blocks[block_index]
        win_start = max(0, start - half)
        win_end = min(proc_rec.get_num_samples(), end + half)
        t_ms = (np.arange(win_start, win_end) - start) * 1000.0 / fs
        span_ms = (end - start) * 1000.0 / fs

        proc = get_traces(proc_rec, win_start, win_end, channel_ids=chosen_ids)
        n_bad = int(np.count_nonzero(proc[(t_ms >= 0) & (t_ms < span_ms)]))

        raw = None
        if raw_rec is not None:
            try:
                raw = get_traces(raw_rec, win_start, win_end, scaled=True, channel_ids=chosen_ids)
            except Exception as err:
                logger.warning(f"Could not read raw traces for the figure: {err}")

        ax_raw = axes[row, 0]
        ax_proc = axes[row, 1]
        if raw is not None:
            for k in range(raw.shape[1]):
                ax_raw.plot(t_ms, raw[:, k] - np.median(raw[:, k]), lw=0.8, alpha=0.8)
        ax_raw.set_ylabel("raw (uV)")
        for k in range(proc.shape[1]):
            ax_proc.plot(t_ms, proc[:, k], lw=0.8, alpha=0.8)
        ax_proc.set_ylabel("saved preprocessed")

        for ax in (ax_raw, ax_proc):
            ax.axvspan(0, span_ms, color="0.85", zorder=0)
            ax.set_xlim(t_ms[0], t_ms[-1])
        ax_raw.set_title(f"block {block_index}: starts {start / fs:.3f} s, {span_ms:.2f} ms", loc="left", fontsize=9)
        ax_proc.set_title(f"nonzero values inside gray span: {n_bad}", loc="left", fontsize=9)
        if row == n_rows - 1:
            ax_raw.set_xlabel("ms from artifact start")
            ax_proc.set_xlabel("ms from artifact start")

    fig.tight_layout()
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


# ----------------------------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------------------------
def main():
    here = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description="Verify that detected artifacts are zeroed in the saved AP data.")
    parser.add_argument("--params", required=True)
    parser.add_argument("--pipeline_script", default=os.path.join(here, PIPELINE_FILENAME))
    parser.add_argument("--probe", type=int, default=None)
    parser.add_argument("--shank", type=int, default=None,
                        help="Shank to check (default: first processed shank if the probe is split by shank)")
    parser.add_argument("--artifacts_csv", default=None,
                        help="Default: <OUTPUT_DIR>/artifact_reward_analysis/imec<N>_artifacts.csv")
    parser.add_argument("--preprocessed_folder", default=None, help="Override the *_si_preprocess folder")
    parser.add_argument("--out_dir", default=None)
    parser.add_argument("--pad_ms", type=float, default=2.0, help="Edge window for CHECK 2")
    parser.add_argument("--flag_sigma", type=float, default=40.0,
                        help="CHECK 2 warns if the edge window exceeds this many noise SDs")
    parser.add_argument("--n_examples", type=int, default=6)
    parser.add_argument("--n_channels", type=int, default=6)
    parser.add_argument("--window_ms", type=float, default=5.0, help="Context shown each side in the figure")
    parser.add_argument("--full_scan", action="store_true", help="CHECK 3: scan the whole file (slow)")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

    pipeline = load_pipeline_module(args.pipeline_script)
    params = pipeline.parse_params_file(args.params)
    probe_idx = args.probe if args.probe is not None else params["PROBE_INDICES"][0]

    shank_idx = args.shank
    if shank_idx is None and probe_idx in params["PROCESS_BY_SHANK_FOR"]:
        if params["SHANKS_TO_PROCESS"] is not None and len(params["SHANKS_TO_PROCESS"]) > 0:
            shank_idx = sorted(params["SHANKS_TO_PROCESS"])[0]
        else:
            shank_idx = 0
        logger.info(f"Probe {probe_idx} is processed by shank; checking shank {shank_idx} (use --shank to change)")

    out_dir = args.out_dir or os.path.join(params["OUTPUT_DIR"], "artifact_reward_analysis")
    os.makedirs(out_dir, exist_ok=True)
    csv_path = args.artifacts_csv or os.path.join(out_dir, f"imec{probe_idx}_artifacts.csv")
    artifacts = pd.read_csv(csv_path)
    blocks = merge_blocks(artifacts["start_sample"].to_numpy(), artifacts["end_sample"].to_numpy())
    logger.info(f"{len(artifacts)} artifact periods -> {len(blocks)} merged zeroed blocks")

    proc_rec, source = open_processed(pipeline, params, probe_idx, shank_idx, blocks, args.preprocessed_folder)
    fs = proc_rec.get_sampling_frequency()
    logger.info(f"Checking: {source}")

    raw_rec = None
    try:
        raw_rec = open_raw(params, probe_idx, shank_idx)
    except Exception as err:
        logger.warning(f"Could not open raw data for the figure: {err}")

    sigma = estimate_noise(proc_rec, blocks, fs)
    table = check_blocks(proc_rec, blocks, fs, args.pad_ms, sigma)

    failed = table[table["n_nonzero_inside"] > 0]
    edge_worst = table[["edge_before_in_sigma", "edge_after_in_sigma"]].max(axis=1)
    warned = table[edge_worst > args.flag_sigma]

    tag = f"imec{probe_idx}" + (f"_shank{shank_idx}" if shank_idx is not None else "")
    lines = []
    lines.append("=" * 78)
    lines.append(f"ZEROING CHECK | {tag} | data checked: {source}")
    lines.append("=" * 78)
    lines.append(f"Merged zeroed blocks checked: {len(table)}   "
                 f"({int(table['n_values_checked'].sum())} channel-samples)")
    lines.append("")
    lines.append("CHECK 1  Are all samples inside the blocks exactly zero?")
    if len(failed) == 0:
        lines.append(f"  PASS: all {len(table)} blocks are exactly zero on every channel.")
    else:
        lines.append(f"  FAIL: {len(failed)} / {len(table)} blocks contain nonzero values:")
        for _, row in failed.head(10).iterrows():
            lines.append(f"    block {int(row['block'])} at {row['start_s']:.4f} s: "
                         f"{int(row['n_nonzero_inside'])} nonzero values, max |x| = {row['max_abs_inside']:.3g}")
    lines.append("")
    lines.append(f"CHECK 2  Edge transients ({args.pad_ms:g} ms either side, noise SD = {sigma:.3g} stored units)")
    lines.append(f"  median edge max: before {np.nanmedian(table['edge_before_in_sigma']):.1f} SD, "
                 f"after {np.nanmedian(table['edge_after_in_sigma']):.1f} SD")
    lines.append(f"  worst edge max:  before {np.nanmax(table['edge_before_in_sigma']):.1f} SD, "
                 f"after {np.nanmax(table['edge_after_in_sigma']):.1f} SD")
    if len(warned) == 0:
        lines.append(f"  OK: no block has an edge window above {args.flag_sigma:g} SD.")
    else:
        lines.append(f"  WARNING: {len(warned)} / {len(table)} blocks have an edge above {args.flag_sigma:g} SD "
                     "(possible leftover artifact just outside the zeroed span, or a real large spike):")
        for _, row in warned.head(10).iterrows():
            lines.append(f"    block {int(row['block'])} at {row['start_s']:.4f} s: "
                         f"before {row['edge_before_in_sigma']:.1f} SD, after {row['edge_after_in_sigma']:.1f} SD")

    if args.full_scan:
        lines.append("")
        lines.append("CHECK 3  Whole-file scan (all-zero samples vs expected)")
        scan = full_scan(proc_rec, blocks, fs)
        lines.append(f"  expected zero samples: {scan['expected_zero_samples']}   "
                     f"observed all-zero samples: {scan['observed_allzero_samples']}")
        lines.append(f"  expected but NOT zero: {scan['expected_but_not_zero']}   "
                     f"zero but NOT expected: {scan['zero_but_not_expected']}")
        if scan["expected_but_not_zero"] == 0 and scan["zero_but_not_expected"] == 0:
            lines.append("  PASS: the zeroed samples are exactly the detected artifacts, nothing more or less.")
        else:
            lines.append("  FAIL/CHECK: mismatch. First unexpected zero times (s): "
                         + ", ".join(f"{t:.4f}" for t in scan["unexpected_times_s"]))
    else:
        lines.append("")
        lines.append("CHECK 3  skipped (add --full_scan to scan the entire file; slow for a full session)")
    lines.append("=" * 78)

    print("\n".join(lines))

    table.to_csv(os.path.join(out_dir, f"{tag}_zeroing_check.csv"), index=False)
    with open(os.path.join(out_dir, f"{tag}_zeroing_check.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")

    examples = pick_examples(blocks, min(args.n_examples, len(blocks)))
    plot_examples(proc_rec, raw_rec, blocks, fs, examples, args.window_ms, args.n_channels,
                  os.path.join(out_dir, f"{tag}_zeroing_examples.png"))
    logger.info(f"Saved results to {out_dir}")


if __name__ == "__main__":
    main()
