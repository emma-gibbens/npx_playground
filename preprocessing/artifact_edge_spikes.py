#!/usr/bin/env python
"""
artifact_edge_spikes.py

Does the zeroing hurt the Kilosort4 sort? This checks your EXISTING sort for the ways zeroing could do
damage, using the Kilosort output (phy format) and the artifact list:

  1. Spikes lost:    how many spikes you would expect to fall inside the zeroed spans (the real yield cost),
                     and how many Kilosort actually placed there.
  2. Edge artifacts: are spikes created right at the edges of the zeroed gaps? (a step from signal to 0, or
                     filter ringing left outside the span, can look like spikes). Shown for all spikes
                     (histograms aligned to gap start and gap end) and per cluster, with a Poisson test of
                     "spikes within --edge_ms of a gap edge" against what that cluster's firing rate predicts.
  3. Bad clusters:   clusters flagged as edge-locked are listed; those are the ones to inspect or drop.

Usage
-----
    python artifact_edge_spikes.py ^
        --ks_folder F:\\KS4_outputs\\M671_2026_08_27_imec0_shank0 ^
        --artifacts_csv D:\\Task4\\M671\\preprocessed\\M671_2026_08_27\\artifact_reward_analysis\\imec0_artifacts.csv

--ks_folder may be the sorter folder or its sorter_output subfolder (the one holding spike_times.npy).
Spike times and artifact times are both in samples of the same recording, so they line up directly.
"""
import argparse
import json
import logging
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

logger = logging.getLogger("artifact_edge_spikes")


def merge_blocks(starts, ends):
    """Merge overlapping/touching [start, end) intervals. Returns two sorted arrays."""
    order = np.argsort(starts)
    block_starts = []
    block_ends = []
    for i in order:
        s = int(starts[i])
        e = int(ends[i])
        if len(block_starts) > 0 and s <= block_ends[-1]:
            if e > block_ends[-1]:
                block_ends[-1] = e
        else:
            block_starts.append(s)
            block_ends.append(e)
    return np.array(block_starts, dtype=np.int64), np.array(block_ends, dtype=np.int64)


def find_sorter_output(folder):
    if os.path.isfile(os.path.join(folder, "spike_times.npy")):
        return folder
    candidate = os.path.join(folder, "sorter_output")
    if os.path.isfile(os.path.join(candidate, "spike_times.npy")):
        return candidate
    raise FileNotFoundError(f"spike_times.npy not found in {folder} or {candidate}")


def poisson_upper_tail(k, mu):
    """P(X >= k) for X ~ Poisson(mu)."""
    try:
        from scipy.stats import poisson
        return float(poisson.sf(k - 1, mu))
    except Exception:
        sd = np.sqrt(max(mu, 1e-12))
        z = (k - 0.5 - mu) / sd
        return float(0.5 * np.exp(-0.717 * z - 0.416 * z * z)) if z > 0 else 1.0


def peri_histogram(spike_times, anchors, edges_samples):
    """Counts of spikes at (spike - anchor) falling in each bin, summed over anchors."""
    counts = np.zeros(len(edges_samples) - 1, dtype=np.int64)
    lo = edges_samples[0]
    hi = edges_samples[-1]
    first = np.searchsorted(spike_times, anchors + lo, side="left")
    last = np.searchsorted(spike_times, anchors + hi, side="left")
    for k in range(len(anchors)):
        if last[k] > first[k]:
            relative = spike_times[first[k]:last[k]] - anchors[k]
            counts += np.histogram(relative, bins=edges_samples)[0]
    return counts


def main():
    parser = argparse.ArgumentParser(description="Check a Kilosort4 sort for zeroing-related damage.")
    parser.add_argument("--ks_folder", required=True)
    parser.add_argument("--artifacts_csv", required=True)
    parser.add_argument("--fs", type=float, default=None, help="Sampling rate (default: from artifacts meta json, else 30000)")
    parser.add_argument("--n_samples", type=int, default=None, help="Recording length in samples (default: meta json)")
    parser.add_argument("--window_ms", type=float, default=10.0, help="Half-width of the peri-gap histograms")
    parser.add_argument("--bin_ms", type=float, default=0.5)
    parser.add_argument("--edge_ms", type=float, default=2.0, help="Edge zone width used for the per-cluster test")
    parser.add_argument("--out_dir", default=None)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

    ks_dir = find_sorter_output(args.ks_folder)
    out_dir = args.out_dir or os.path.dirname(os.path.abspath(args.artifacts_csv))
    os.makedirs(out_dir, exist_ok=True)

    fs = args.fs
    n_samples = args.n_samples
    meta_guess = args.artifacts_csv.replace("_artifacts.csv", "_artifacts_meta.json")
    if os.path.isfile(meta_guess):
        with open(meta_guess) as f:
            meta = json.load(f)
        if fs is None:
            fs = float(meta["fs"])
        if n_samples is None:
            n_samples = int(meta["n_samples"])
    if fs is None:
        fs = 30000.0

    spike_times = np.load(os.path.join(ks_dir, "spike_times.npy")).astype(np.int64).ravel()
    spike_clusters = np.load(os.path.join(ks_dir, "spike_clusters.npy")).astype(np.int64).ravel()
    order = np.argsort(spike_times, kind="stable")
    spike_times = spike_times[order]
    spike_clusters = spike_clusters[order]
    if n_samples is None:
        n_samples = int(spike_times.max()) + 1
    total_s = n_samples / fs

    labels = {}
    label_file = os.path.join(ks_dir, "cluster_KSLabel.tsv")
    if os.path.isfile(label_file):
        table = pd.read_csv(label_file, sep="\t")
        for _, row in table.iterrows():
            labels[int(row["cluster_id"])] = str(row["KSLabel"])

    artifacts = pd.read_csv(args.artifacts_csv)
    block_starts, block_ends = merge_blocks(artifacts["start_sample"].to_numpy(), artifacts["end_sample"].to_numpy())
    n_blocks = len(block_starts)
    zeroed_samples = int((block_ends - block_starts).sum())
    zeroed_s = zeroed_samples / fs
    n_spikes = len(spike_times)
    logger.info(f"{n_spikes} spikes, {len(np.unique(spike_clusters))} clusters, {n_blocks} zeroed blocks, "
                f"{zeroed_s * 1000:.1f} ms zeroed of {total_s:.1f} s")

    lines = []
    lines.append("=" * 78)
    lines.append(f"ZEROING vs SORT | {ks_dir}")
    lines.append("=" * 78)

    # ---- 1. spikes lost ------------------------------------------------------------------------
    idx = np.searchsorted(block_starts, spike_times, side="right") - 1
    inside = np.zeros(n_spikes, dtype=bool)
    valid = idx >= 0
    inside[valid] = spike_times[valid] < block_ends[idx[valid]]
    expected_lost = n_spikes / total_s * zeroed_s / max(1e-12, (1 - zeroed_s / total_s))
    lines.append("1. SPIKES LOST TO ZEROING")
    lines.append(f"   zeroed time: {zeroed_s * 1000:.1f} ms in {n_blocks} blocks = {100 * zeroed_s / total_s:.4f} % of the recording")
    lines.append(f"   spikes expected inside zeroed spans if firing were uniform: {expected_lost:.1f} "
                 f"({100 * expected_lost / (n_spikes + expected_lost):.4f} % of all spikes)")
    lines.append(f"   spikes Kilosort actually placed inside zeroed spans: {int(inside.sum())}")
    lines.append("")

    # ---- 2. edge spikes (all clusters pooled) --------------------------------------------------
    edge = int(round(args.edge_ms * fs / 1000.0))
    next_i = np.searchsorted(block_starts, spike_times, side="right")
    before_start = np.zeros(n_spikes, dtype=bool)
    has_next = next_i < n_blocks
    gap = np.zeros(n_spikes, dtype=np.int64)
    gap[has_next] = block_starts[next_i[has_next]] - spike_times[has_next]
    before_start[has_next] = (gap[has_next] > 0) & (gap[has_next] <= edge)

    prev_i = np.searchsorted(block_ends, spike_times, side="right") - 1
    after_end = np.zeros(n_spikes, dtype=bool)
    has_prev = prev_i >= 0
    since = np.zeros(n_spikes, dtype=np.int64)
    since[has_prev] = spike_times[has_prev] - block_ends[prev_i[has_prev]]
    after_end[has_prev] = (since[has_prev] >= 0) & (since[has_prev] < edge)

    edge_spike = (before_start | after_end) & ~inside
    edge_total_s = n_blocks * 2 * args.edge_ms / 1000.0
    rate_all = n_spikes / total_s
    expected_edge = rate_all * edge_total_s
    n_edge = int(edge_spike.sum())
    lines.append(f"2. SPIKES AT GAP EDGES (within {args.edge_ms:g} ms before a gap or after a gap)")
    lines.append(f"   observed {n_edge}   expected from the overall rate {expected_edge:.1f}   "
                 f"ratio {n_edge / max(expected_edge, 1e-9):.2f}  (p = {poisson_upper_tail(n_edge, expected_edge):.3g})")

    # ---- 3. per cluster ------------------------------------------------------------------------
    cluster_ids, inverse = np.unique(spike_clusters, return_inverse=True)
    n_per_cluster = np.bincount(inverse, minlength=len(cluster_ids))
    n_edge_per_cluster = np.bincount(inverse[edge_spike], minlength=len(cluster_ids))
    rows = []
    for k in range(len(cluster_ids)):
        mu = n_per_cluster[k] / total_s * edge_total_s
        k_obs = int(n_edge_per_cluster[k])
        p = poisson_upper_tail(k_obs, mu) if k_obs > 0 else 1.0
        rows.append(dict(cluster_id=int(cluster_ids[k]), ks_label=labels.get(int(cluster_ids[k]), ""),
                         n_spikes=int(n_per_cluster[k]), n_edge_spikes=k_obs, expected_edge=mu,
                         ratio=k_obs / max(mu, 1e-9), p_upper=p))
    per_cluster = pd.DataFrame(rows)
    n_tests = len(per_cluster)
    per_cluster["p_bonferroni"] = np.minimum(1.0, per_cluster["p_upper"] * n_tests)
    per_cluster["edge_locked"] = (per_cluster["p_bonferroni"] < 0.05) & (per_cluster["n_edge_spikes"] >= 3)
    flagged = per_cluster[per_cluster["edge_locked"]].sort_values("ratio", ascending=False)
    lines.append("")
    lines.append("3. PER-CLUSTER EDGE LOCKING (Poisson test vs each cluster's own rate, Bonferroni corrected)")
    lines.append(f"   clusters tested: {n_tests};  flagged edge-locked: {len(flagged)}")
    if len(flagged) > 0:
        lines.append("   cluster  label  n_spikes  edge_spikes  expected  ratio")
        for _, row in flagged.head(15).iterrows():
            lines.append(f"   {int(row['cluster_id']):7d}  {row['ks_label']:>5s}  {int(row['n_spikes']):8d}  "
                         f"{int(row['n_edge_spikes']):11d}  {row['expected_edge']:8.2f}  {row['ratio']:6.1f}")
        share = flagged["n_edge_spikes"].sum() / max(1, n_edge)
        lines.append(f"   flagged clusters account for {100 * share:.0f} % of all edge spikes")
    lines.append("=" * 78)
    print("\n".join(lines))

    # ---- figure --------------------------------------------------------------------------------
    half = int(round(args.window_ms * fs / 1000.0))
    bin_samples = max(1, int(round(args.bin_ms * fs / 1000.0)))
    edges = np.arange(-half, half + bin_samples, bin_samples)
    centers_ms = 0.5 * (edges[:-1] + edges[1:]) * 1000.0 / fs
    hist_start = peri_histogram(spike_times, block_starts, edges)
    hist_end = peri_histogram(spike_times, block_ends, edges)

    fig, axes = plt.subplots(1, 3, figsize=(17, 4.4))
    for ax, hist, title, span in ((axes[0], hist_start, "A. Spikes aligned to gap START", (0, None)),
                                  (axes[1], hist_end, "B. Spikes aligned to gap END", (None, 0))):
        ax.bar(centers_ms, hist, width=args.bin_ms * 0.95, color="0.4")
        far = np.abs(centers_ms) >= args.window_ms * 0.5
        far_values = hist[far & (hist > -1)]
        if len(far_values) > 0:
            ax.axhline(far_values.mean(), color="crimson", ls="--", lw=1, label="mean of far bins")
        ax.axvline(0, color="tab:blue", lw=1)
        ax.set_xlabel("ms from gap " + ("start" if span[0] == 0 else "end"))
        ax.set_ylabel("Spikes (all clusters)")
        ax.set_title(title + "  (zeroed side is " + ("right" if span[0] == 0 else "left") + " of 0)",
                     loc="left", fontsize=9)
        ax.legend(frameon=False)

    ax = axes[2]
    ax.scatter(per_cluster["n_spikes"], per_cluster["ratio"], s=14, c="0.5", label="clusters")
    if len(flagged) > 0:
        ax.scatter(flagged["n_spikes"], flagged["ratio"], s=26, c="crimson", label="edge-locked")
    ax.axhline(1.0, color="k", lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel("Spikes in cluster")
    ax.set_ylabel("Edge spikes / expected")
    ax.set_title("C. Per-cluster edge locking", loc="left", fontsize=9)
    ax.legend(frameon=False)
    fig.tight_layout()

    tag = os.path.basename(os.path.normpath(args.ks_folder if os.path.basename(args.ks_folder) != "sorter_output"
                                            else os.path.dirname(args.ks_folder)))
    fig.savefig(os.path.join(out_dir, f"{tag}_edge_spikes.png"), dpi=170)
    plt.close(fig)
    per_cluster.to_csv(os.path.join(out_dir, f"{tag}_edge_spikes_per_cluster.csv"), index=False)
    with open(os.path.join(out_dir, f"{tag}_edge_spikes.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    logger.info(f"Saved results to {out_dir}")


if __name__ == "__main__":
    main()
