#!/usr/bin/env python
"""
artifacts_vs_reward.py

Detect saturation artifacts with the SAME code/parameters as extract_lfp_and_run_ks4_update.py,
then characterise them (how much recording gets zeroed) and relate them to reward timing.

What it does
------------
1. Reads the params .txt file used by the pipeline (so CATGT_OUTPUT_DIR, ARTIFACT_PROPORTION,
   ARTIFACT_SATURATION_THRESHOLD_UV, SI_NUM_JOBS, CHUNK_DURATION are identical to the real run).
2. Loads imec{probe}.ap from CATGT_OUTPUT_DIR with si.read_spikeglx and calls the pipeline's own
   detect_artifacts(), i.e. si.detect_artifact_periods(method='saturation'). The periods it returns
   are exactly what si.silence_periods(..., mode='zeros') zeroes out in the pipeline.
   (Detection is slow on a full session, so the result is cached; use --reuse to skip it next time.)
3. Reports artifact length (mean / median / max / min), total zeroed time and % of recording.
   Reported for the raw detected periods AND for merged periods (overlapping/touching periods
   are one continuous zeroed block, which is what actually ends up in your data).
4. Relates artifacts to reward times (XD1_ON.txt, seconds):
     - reward-triggered histogram of artifact onsets
     - fraction of trials in which the data at each peri-reward time is zeroed
     - per-artifact latency from the previous reward
     - reward-locked artifacts: onset within [-lock_pre, +lock_post] (default -100 to +500 ms) of a reward,
       with their count, onset timing, length, zeroed time and a permutation test
     - circular-shift permutation test: are artifacts more/less common after reward than chance?

Usage (from an environment that has spikeinterface installed)
-------------------------------------------------------------
    python artifacts_vs_reward.py --params D:\\path\\to\\M671_2026_08_27_lfp_ks4_params.txt

    # re-plot / re-analyse without redoing detection
    python artifacts_vs_reward.py --params ...txt --reuse

    # custom windows
    python artifacts_vs_reward.py --params ...txt --reuse --pre 2 --post 5 --test_window 0 1.0

Outputs (default: <OUTPUT_DIR>\\artifact_reward_analysis)
-------------------------------------------------------
    imec0_artifacts.csv          every artifact: start/end (samples, s), duration, reward relationships
    imec0_artifacts_meta.json    sampling rate / recording length (used by --reuse)
    imec0_artifact_summary.txt   the numbers printed at the end
    imec0_artifacts_vs_reward.png / .pdf
"""
import argparse
import importlib.util
import json
import logging
import os
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

PIPELINE_FILENAME = "extract_lfp_and_run_ks4_update.py"
logger = logging.getLogger("artifacts_vs_reward")


# ----------------------------------------------------------------------------------------------
# Loading the pipeline + detecting artifacts
# ----------------------------------------------------------------------------------------------
def load_pipeline_module(path):
    """Import the existing pipeline script as a module (its main() is guarded by __name__)."""
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"Could not find the pipeline script at {path}. Pass --pipeline_script explicitly."
        )
    spec = importlib.util.spec_from_file_location("ks4_pipeline", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def detect_with_pipeline(pipeline, params, probe_idx):
    """Run the pipeline's own artifact detection on the full (unsplit) AP stream."""
    import spikeinterface.full as si

    p = dict(params)
    if not (p["ARTIFACT_REJECTION_AP"] or p["ARTIFACT_REJECTION_LFP"]):
        logger.warning("Artifact rejection is OFF in the params file; forcing detection for this analysis.")
        p["ARTIFACT_REJECTION_AP"] = True

    logger.info(f"Reading imec{probe_idx}.ap from {p['CATGT_OUTPUT_DIR']} ...")
    rec = si.read_spikeglx(p["CATGT_OUTPUT_DIR"], stream_name=f"imec{probe_idx}.ap")
    logger.info(
        f"Recording: {rec.get_num_channels()} channels, fs={rec.get_sampling_frequency():.3f} Hz, "
        f"{rec.get_num_samples(0) / rec.get_sampling_frequency():.1f} s, {rec.get_num_segments()} segment(s)"
    )
    # same function, same params as the real run (proportion, saturation threshold, jobs, chunk size)
    artifacts = pipeline.detect_artifacts(rec, probe_idx, logger, p)
    return rec, artifacts


def artifacts_to_samples(artifacts):
    """Convert SpikeInterface's artifact-period structured array into (start, end) sample arrays."""
    arr = np.asarray(artifacts)
    if arr.size == 0:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)
    names = arr.dtype.names
    if names is None:
        raise ValueError(f"Unexpected artifact format (dtype={arr.dtype}); expected a structured array.")

    def pick(*candidates, contains=None):
        for c in candidates:
            if c in names:
                return c
        if contains:
            for n in names:
                if contains in n.lower():
                    return n
        raise ValueError(f"Could not find a '{candidates[0]}' field in artifact dtype names {names}")

    start_f = pick("start_sample_index", "start_sample", contains="start")
    end_f = pick("end_sample_index", "end_sample", contains="end")

    if "segment_index" in names:
        keep = arr["segment_index"] == 0
        if not np.all(keep):
            logger.warning("Recording has more than one segment; only segment 0 is analysed.")
        arr = arr[keep]
    return arr[start_f].astype(np.int64), arr[end_f].astype(np.int64)


# ----------------------------------------------------------------------------------------------
# Interval helpers
# ----------------------------------------------------------------------------------------------
def merge_intervals(starts, ends):
    """Merge overlapping / touching [start, end) intervals. Returns (merged_starts, merged_ends)."""
    if len(starts) == 0:
        return np.array([]), np.array([])
    order = np.argsort(starts, kind="stable")
    s, e = starts[order], ends[order]
    running_end = np.maximum.accumulate(e)
    new_block = np.r_[0, np.where(s[1:] > running_end[:-1])[0] + 1]
    return s[new_block], np.maximum.reduceat(e, new_block)


def inside_intervals(times, merged_s, merged_e):
    """Boolean array: is each time inside any merged interval?"""
    if len(merged_s) == 0:
        return np.zeros(np.shape(times), dtype=bool)
    idx = np.searchsorted(merged_s, times, side="right") - 1
    ok = idx >= 0
    out = np.zeros(np.shape(times), dtype=bool)
    out[ok] = times[ok] < merged_e[idx[ok]]
    return out


def zeroed_overlap(win_start, win_end, merged_s, merged_e):
    """Seconds of each window [win_start, win_end) that fall inside merged intervals."""
    out = np.zeros(len(win_start))
    for i, (a, b) in enumerate(zip(win_start, win_end)):
        lo = max(np.searchsorted(merged_e, a, side="right"), 0)
        hi = np.searchsorted(merged_s, b, side="left")
        for k in range(lo, hi):
            out[i] += max(0.0, min(b, merged_e[k]) - max(a, merged_s[k]))
    return out


# ----------------------------------------------------------------------------------------------
# Peri-reward statistics
# ----------------------------------------------------------------------------------------------
def peri_counts(onsets_sorted, rewards, edges, tw):
    """
    onsets_sorted : sorted artifact onset times (s)
    rewards       : reward times (s)
    edges         : histogram bin edges relative to reward
    tw            : (t0, t1) test window relative to reward
    Returns (counts per bin summed over rewards, number of onsets in the test window for each reward).
    """
    idx = np.searchsorted(onsets_sorted, rewards[:, None] + edges[None, :], side="left")
    counts = np.diff(idx, axis=1).sum(axis=0)
    lo = np.searchsorted(onsets_sorted, rewards + tw[0], side="left")
    hi = np.searchsorted(onsets_sorted, rewards + tw[1], side="left")
    return counts, hi - lo


def permutation_test(onsets, rewards, total_dur, edges_list, tw, n_shuffles, min_shift, seed=0):
    """Circularly shift ALL artifact onsets by a random amount (keeps artifact clustering and reward
    spacing intact) and recompute the peri-reward statistics (one histogram per edge set in edges_list)."""
    rng = np.random.default_rng(seed)
    null_hists = [np.zeros((n_shuffles, len(e) - 1)) for e in edges_list]
    null_frac = np.zeros(n_shuffles)   # fraction of rewards with >=1 artifact onset in test window
    null_total = np.zeros(n_shuffles)  # total artifact onsets in test windows
    hi_shift = total_dur - min_shift
    if hi_shift <= min_shift:
        raise ValueError("Recording too short for the requested minimum shift.")
    for k in range(n_shuffles):
        shift = rng.uniform(min_shift, hi_shift)
        sh = np.sort((onsets + shift) % total_dur)
        for j, edges in enumerate(edges_list):
            h, per_reward = peri_counts(sh, rewards, edges, tw)
            null_hists[j][k] = h
            if j == 0:
                null_frac[k] = np.mean(per_reward > 0)
                null_total[k] = per_reward.sum()
    return null_hists, null_frac, null_total


def perm_p(obs, null):
    """Two-sided permutation p-value (+1 correction) and direction."""
    ge = (1 + np.sum(null >= obs)) / (len(null) + 1)
    le = (1 + np.sum(null <= obs)) / (len(null) + 1)
    return min(1.0, 2 * min(ge, le)), ("more than chance" if ge < le else "less than chance")


def zeroed_bin_hits(rewards, edges, merged_s, merged_e):
    """
    (n_rewards x n_bins) boolean: does ANY zeroed span overlap bin [edge_i, edge_i+1) after each reward?
    (Checking only the bin centre would miss ~1 ms artifacts almost every time.)
    merged_s / merged_e must be the non-overlapping, sorted output of merge_intervals().
    """
    a = rewards[:, None] + edges[None, :-1]
    b = rewards[:, None] + edges[None, 1:]
    n_overlap = np.searchsorted(merged_s, b, side="left") - np.searchsorted(merged_e, a, side="right")
    return n_overlap > 0


# ----------------------------------------------------------------------------------------------
# Main analysis
# ----------------------------------------------------------------------------------------------
def describe_durations(d_s):
    if len(d_s) == 0:
        return {k: np.nan for k in ("mean", "median", "max", "min", "std")}
    return dict(mean=np.mean(d_s), median=np.median(d_s), max=np.max(d_s), min=np.min(d_s), std=np.std(d_s))


def fmt_ms(x):
    return "n/a" if not np.isfinite(x) else f"{x * 1000:.2f} ms"


def make_edges(lo, hi, width):
    return np.linspace(lo, hi, int(round((hi - lo) / width)) + 1)


def analyse(start_s, end_s, rewards, total_dur, args, probe_idx):
    lines = []
    P = lines.append

    dur = end_s - start_s
    n_art = len(start_s)
    ms, me = merge_intervals(start_s, end_s)
    merged_dur = me - ms
    total_zeroed = merged_dur.sum()

    P("=" * 78)
    P(f"ARTIFACTS vs REWARD  |  probe {probe_idx}  |  recording length {total_dur:.1f} s ({total_dur / 60:.1f} min)")
    P("=" * 78)
    P(f"Detected artifact periods: {n_art}   ({n_art / (total_dur / 60):.2f} per minute)")
    P("")
    P("Artifact length (raw detected periods)")
    d = describe_durations(dur)
    P(f"  mean   {fmt_ms(d['mean'])}   (SD {fmt_ms(d['std'])})")
    P(f"  median {fmt_ms(d['median'])}")
    P(f"  max    {fmt_ms(d['max'])}" + (f"   (starts at {start_s[np.argmax(dur)]:.3f} s)" if n_art else ""))
    P(f"  min    {fmt_ms(d['min'])}")
    P("")
    P("Merged zeroed blocks (overlapping/touching periods combined = what is really zeroed)")
    dm = describe_durations(merged_dur)
    P(f"  blocks {len(ms)}")
    P(f"  mean   {fmt_ms(dm['mean'])}   median {fmt_ms(dm['median'])}   max {fmt_ms(dm['max'])}")
    P(f"  TOTAL zeroed: {total_zeroed:.3f} s  =  {100 * total_zeroed / total_dur:.4f} % of the recording")
    P("")

    result = dict(merged_s=ms, merged_e=me)
    if len(rewards) == 0 or n_art == 0:
        P("No rewards or no artifacts -> skipping reward analysis.")
        return lines, result

    # --- reward relationships -----------------------------------------------------------------
    tw = tuple(args.test_window)
    edges_w = make_edges(-args.pre, args.post, args.bin)                   # wide view
    edges_z = make_edges(-args.zoom_pre, args.zoom_post, args.zoom_bin)    # zoomed view (ms scale)
    centers_w = 0.5 * (edges_w[:-1] + edges_w[1:])
    centers_z = 0.5 * (edges_z[:-1] + edges_z[1:])

    rew_in_art = inside_intervals(rewards, ms, me)

    # per-artifact latency to neighbouring rewards
    prev_idx = np.searchsorted(rewards, start_s, side="right") - 1
    t_since_prev = np.where(prev_idx >= 0, start_s - rewards[np.clip(prev_idx, 0, None)], np.nan)
    next_idx = np.searchsorted(rewards, start_s, side="left")
    t_to_next = np.where(next_idx < len(rewards), rewards[np.clip(next_idx, None, len(rewards) - 1)] - start_s, np.nan)
    in_window = (t_since_prev >= tw[0]) & (t_since_prev < tw[1])

    sorted_on = np.sort(start_s)
    obs_hist_w, obs_per_reward = peri_counts(sorted_on, rewards, edges_w, tw)
    obs_hist_z, _ = peri_counts(sorted_on, rewards, edges_z, tw)
    obs_frac = np.mean(obs_per_reward > 0)
    obs_total = obs_per_reward.sum()

    min_shift = max(args.pre + args.post + 5.0, 30.0)
    (null_w, null_z), null_frac, null_total = permutation_test(
        start_s, rewards, total_dur, [edges_w, edges_z], tw, args.n_shuffles, min_shift, seed=args.seed
    )
    p_frac, dir_frac = perm_p(obs_frac, null_frac)
    p_tot, dir_tot = perm_p(obs_total, null_total)

    # zeroed-data coverage around reward (any zeroed sample inside each bin)
    cov_hits = zeroed_bin_hits(rewards, edges_z, ms, me)
    coverage = cov_hits.mean(axis=0)
    win_zeroed = zeroed_overlap(rewards + tw[0], rewards + tw[1], ms, me)
    win_len = tw[1] - tw[0]

    P(f"Reward timing  ({len(rewards)} rewards; first {rewards[0]:.3f} s, last {rewards[-1]:.3f} s)")
    P(f"  Rewards that land INSIDE a zeroed block: {rew_in_art.sum()} / {len(rewards)}")
    P(f"  Test window = {tw[0]:g} to {tw[1]:g} s after reward onset")
    P(f"  Mean fraction of that window zeroed: {100 * win_zeroed.mean() / win_len:.3f} %   "
      f"(whole-recording baseline {100 * total_zeroed / total_dur:.4f} %)")
    P(f"  Max zeroed in any one reward window: {win_zeroed.max() * 1000:.1f} ms")
    P(f"  Rewards with >=1 artifact onset in window: {obs_frac * 100:.1f} %  "
      f"(chance {100 * null_frac.mean():.1f} % +/- {100 * null_frac.std():.1f}; "
      f"permutation p = {p_frac:.4f}, {dir_frac})")
    P(f"  Artifact onsets in window (summed over rewards): {obs_total}  "
      f"(chance {null_total.mean():.1f} +/- {null_total.std():.1f}; p = {p_tot:.4f}, {dir_tot})")
    n_win = int(in_window.sum())
    P(f"  Artifacts whose onset falls in the post-reward window: {n_win} / {n_art} ({100 * n_win / n_art:.1f} %)")
    if n_win and (n_art - n_win):
        P(f"  Mean length  in-window {fmt_ms(dur[in_window].mean())}  vs  other {fmt_ms(dur[~in_window].mean())}")
    valid = np.isfinite(t_since_prev)
    if valid.any():
        P(f"  Latency from previous reward to artifact onset: median {np.nanmedian(t_since_prev):.3f} s "
          f"(25-75%: {np.nanpercentile(t_since_prev, 25):.3f}-{np.nanpercentile(t_since_prev, 75):.3f} s)")

    # ---- reward-locked artifacts: onset within [-lock_pre, +lock_post) of the NEAREST reward ------
    lk0, lk1 = -args.lock_pre, args.lock_post
    since = np.where(np.isfinite(t_since_prev), t_since_prev, np.inf)
    until = np.where(np.isfinite(t_to_next), t_to_next, np.inf)
    near_t = np.where(since <= until, since, -until)      # signed: + = after reward, - = before reward
    locked = (near_t >= lk0) & (near_t < lk1)
    n_lock = int(locked.sum())

    lock_hist, lock_per_reward = peri_counts(sorted_on, rewards, edges_z, (lk0, lk1))
    lock_obs_frac = np.mean(lock_per_reward > 0)
    _, lock_null_frac, lock_null_total = permutation_test(
        start_s, rewards, total_dur, [edges_z], (lk0, lk1), args.n_shuffles, min_shift, seed=args.seed)
    lock_p, lock_dir = perm_p(lock_obs_frac, lock_null_frac)
    lock_zeroed = zeroed_overlap(rewards + lk0, rewards + lk1, ms, me)
    lock_len = lk1 - lk0

    P("")
    P(f"REWARD-LOCKED artifacts: onset from {lk0 * 1000:+.0f} ms to {lk1 * 1000:+.0f} ms around reward onset")
    P(f"  Artifacts in window: {n_lock} / {n_art} ({100 * n_lock / n_art:.1f} % of all artifacts)")
    P(f"  Rewards with >=1 artifact in window: {int((lock_per_reward > 0).sum())} / {len(rewards)} "
      f"({100 * lock_obs_frac:.1f} %)   chance {100 * lock_null_frac.mean():.1f} % +/- {100 * lock_null_frac.std():.1f}; "
      f"permutation p = {lock_p:.4f}, {lock_dir}")
    P(f"  Artifact onsets in window per reward: mean {lock_per_reward.mean():.2f}, max {int(lock_per_reward.max())}   "
      f"(chance total {lock_null_total.mean():.1f} vs observed {int(lock_per_reward.sum())})")
    if n_lock:
        lt = near_t[locked]
        P(f"  Onset relative to reward: median {np.median(lt) * 1000:+.1f} ms, "
          f"IQR {np.percentile(lt, 25) * 1000:+.1f} to {np.percentile(lt, 75) * 1000:+.1f} ms, "
          f"range {lt.min() * 1000:+.1f} to {lt.max() * 1000:+.1f} ms")
        P(f"  Onset timing (cumulative after reward): " + ",  ".join(
            f"<{c * 1000:.0f} ms: {int(((lt >= 0) & (lt < c)).sum())}" for c in (0.05, 0.1, 0.25, 0.5) if c <= lk1 + 1e-9)
          + f",  BEFORE logged reward: {int((lt < 0).sum())}")
        P(f"  Length of reward-locked artifacts: mean {fmt_ms(dur[locked].mean())}, "
          f"median {fmt_ms(np.median(dur[locked]))}, max {fmt_ms(dur[locked].max())}   "
          f"(all other artifacts: mean {fmt_ms(dur[~locked].mean()) if (~locked).any() else 'n/a'})")
    P(f"  Zeroed time inside the {lock_len * 1000:.0f} ms window: mean {lock_zeroed.mean() * 1000:.3f} ms per reward "
      f"({100 * lock_zeroed.mean() / lock_len:.3f} % of window), max {lock_zeroed.max() * 1000:.2f} ms, "
      f"total {lock_zeroed.sum() * 1000:.1f} ms across all rewards")
    P("")
    pk = np.argmax(obs_hist_z)
    P(f"  Peak of fine ({args.zoom_bin * 1000:g} ms) reward-triggered histogram: {centers_z[pk] * 1000:+.1f} ms "
      f"({obs_hist_z[pk]} onsets vs {null_z.mean(axis=0)[pk]:.2f} expected by chance)")
    P("=" * 78)

    result.update(dict(
        edges_w=edges_w, centers_w=centers_w, obs_hist_w=obs_hist_w, null_w=null_w,
        edges_z=edges_z, centers_z=centers_z, obs_hist_z=obs_hist_z, null_z=null_z,
        null_frac=null_frac, obs_frac=obs_frac, p_frac=p_frac, coverage=coverage, tw=tw,
        rew_in_art=rew_in_art, t_since_prev=t_since_prev, t_to_next=t_to_next, in_window=in_window, dur=dur,
        near_t=near_t, locked=locked))
    return lines, result


# ----------------------------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------------------------
def peth_panel(ax, centers, edges, obs_hist, null_hist, n_rew, xscale, xlabel, tw, title):
    bin_w = edges[1] - edges[0]
    rate = obs_hist / (n_rew * bin_w)
    null_rate = null_hist / (n_rew * bin_w)
    lo_b, hi_b = np.percentile(null_rate, [2.5, 97.5], axis=0)
    ax.fill_between(centers * xscale, lo_b, hi_b, color="0.8", step="mid", label="95% of shuffled data")
    ax.step(centers * xscale, rate, where="mid", color="crimson", lw=1.5, label="observed")
    ax.axvline(0, color="tab:blue", lw=1.0, zorder=0)
    ax.set_xlim(edges[0] * xscale, edges[-1] * xscale)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Artifact onsets per s per reward")
    ax.set_title(title, loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8)


def make_figure(start_s, end_s, rewards, total_dur, res, args, probe_idx, summary_lines, out_png, out_pdf):
    dur = end_s - start_s
    fig = plt.figure(figsize=(15, 20))
    gs = GridSpec(5, 2, figure=fig, height_ratios=[0.9, 1.3, 1.3, 1.2, 1.1], hspace=0.42, wspace=0.22)

    # A: whole-session timeline
    ax = fig.add_subplot(gs[0, :])
    ax.vlines(start_s / 60, 0, 0.45, color="crimson", lw=0.8, label=f"artifact onset (n={len(start_s)})")
    ax.vlines(rewards / 60, 0.55, 1.0, color="tab:blue", lw=1.0, label=f"reward (n={len(rewards)})")
    ax.set_xlim(0, total_dur / 60)
    ax.set_ylim(0, 1.3)
    ax.set_yticks([])
    ax.set_xlabel("Time in recording (min)")
    ax.set_title(f"A. Probe {probe_idx}: artifact onsets and reward times across the session", loc="left")
    ax.legend(loc="upper right", ncol=2, frameon=False)

    # B: duration histogram
    ax = fig.add_subplot(gs[1, 0])
    if len(dur):
        dms = dur * 1000
        bins = np.logspace(np.log10(max(dms.min(), 1e-3)), np.log10(dms.max() * 1.05 + 1e-3), 40)
        ax.hist(dms, bins=bins, color="0.35")
        ax.set_xscale("log")
        ax.axvline(dms.mean(), color="tab:orange", lw=2, label=f"mean {dms.mean():.2f} ms")
        ax.axvline(dms.max(), color="crimson", lw=2, ls="--", label=f"max {dms.max():.2f} ms")
    ax.set_xlabel("Artifact length (ms, log scale)")
    ax.set_ylabel("Count")
    ax.set_title("B. Length of detected (zeroed) artifacts", loc="left")
    ax.legend(frameon=False)

    if "centers_z" not in res:  # no reward analysis possible
        fig.savefig(out_png, dpi=200, bbox_inches="tight")
        fig.savefig(out_pdf, bbox_inches="tight")
        plt.close(fig)
        return

    tw = res["tw"]
    n_rew = len(rewards)
    zpre, zpost = args.zoom_pre, args.zoom_post

    # C: zoomed raster aligned to reward (ms scale so ~1 ms artifacts are visible)
    ax = fig.add_subplot(gs[1, 1])
    ax.axvline(0, color="tab:blue", lw=1.0, zorder=0)
    n_with = 0
    for i, r in enumerate(rewards):
        lo = np.searchsorted(end_s, r - zpre, side="left")   # artifacts that end after window start
        hi = np.searchsorted(start_s, r + zpost, side="left")
        if hi > lo:
            n_with += 1
            a = (start_s[lo:hi] - r) * 1000
            b = (end_s[lo:hi] - r) * 1000
            ax.plot(a, np.full(hi - lo, i), "o", color="crimson", ms=4, mew=0, zorder=3)  # onset (always visible)
            ax.hlines(np.full(hi - lo, i), a, b, color="k", lw=2.0, zorder=2)              # zeroed span (to scale)
    ax.set_xlim(-zpre * 1000, zpost * 1000)
    ax.set_ylim(n_rew - 0.5, -0.5)
    ax.set_xlabel("Time from reward (ms)")
    ax.set_ylabel("Reward #")
    ax.set_title(f"C. Artifacts per reward: {n_with}/{n_rew} rewards have one (dot = onset, bar = span)",
                 loc="left", fontsize=10)

    # D: wide PETH
    ax = fig.add_subplot(gs[2, 0])
    peth_panel(ax, res["centers_w"], res["edges_w"], res["obs_hist_w"], res["null_w"], n_rew, 1.0,
               "Time from reward (s)", tw, f"D. Reward-triggered artifact-onset rate, wide view ({args.bin * 1000:g} ms bins)")
    ax.axvspan(tw[0], tw[1], color="tab:blue", alpha=0.08)

    # E: zoomed PETH
    ax = fig.add_subplot(gs[2, 1])
    peth_panel(ax, res["centers_z"], res["edges_z"], res["obs_hist_z"], res["null_z"], n_rew, 1000.0,
               "Time from reward (ms)", tw, f"E. Same, zoomed ({args.zoom_bin * 1000:g} ms bins)")

    # F: coverage (any zeroed sample in each bin)
    ax = fig.add_subplot(gs[3, 0])
    ax.step(res["centers_z"] * 1000, 100 * res["coverage"], where="mid", color="k", lw=1.5)
    ax.axvline(0, color="tab:blue", lw=1.0, zorder=0)
    ax.set_xlim(-zpre * 1000, zpost * 1000)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Time from reward (ms)")
    ax.set_ylabel("% of rewards")
    ax.set_title(f"F. % of rewards with ANY zeroed samples in each {args.zoom_bin * 1000:g} ms bin", loc="left", fontsize=10)

    # G: permutation null
    ax = fig.add_subplot(gs[3, 1])
    ax.hist(100 * res["null_frac"], bins=30, color="0.7", label="shuffled")
    ax.axvline(100 * res["obs_frac"], color="crimson", lw=2,
               label=f"observed {100 * res['obs_frac']:.1f}%  (p={res['p_frac']:.3g})")
    ax.set_xlabel(f"% of rewards with an artifact onset {tw[0]:g} to {tw[1]:g} s after reward")
    ax.set_ylabel("Shuffles")
    ax.set_title("G. Permutation test (artifacts circularly shifted in time)", loc="left", fontsize=10)
    ax.legend(frameon=False)

    # H: summary text
    ax = fig.add_subplot(gs[4, :])
    ax.axis("off")
    ax.text(0, 1, "\n".join(l for l in summary_lines if not set(l) <= {"="}), va="top", ha="left",
            family="monospace", fontsize=8, transform=ax.transAxes)

    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------------------------
def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description="Relate zeroed artifacts to reward timing (SpikeInterface).")
    ap.add_argument("--params", required=True, help="The same params .txt used for extract_lfp_and_run_ks4_update.py")
    ap.add_argument("--pipeline_script", default=os.path.join(here, PIPELINE_FILENAME),
                    help="Path to extract_lfp_and_run_ks4_update.py (default: same folder as this script)")
    ap.add_argument("--probe", type=int, default=None, help="Probe index (default: first in PROBE_INDICES)")
    ap.add_argument("--reward_file", default=None, help="Reward times in s (default: <OUTPUT_DIR>/XD1_ON.txt)")
    ap.add_argument("--reward_offset", type=float, default=0.0,
                    help="Seconds added to reward times if they are not on the same clock as the .ap file")
    ap.add_argument("--out_dir", default=None, help="Default: <OUTPUT_DIR>/artifact_reward_analysis")
    ap.add_argument("--reuse", action="store_true", help="Reuse cached detection instead of re-running it")
    ap.add_argument("--pre", type=float, default=2.0, help="Seconds before reward shown in peri-reward plots")
    ap.add_argument("--post", type=float, default=5.0, help="Seconds after reward shown in peri-reward plots")
    ap.add_argument("--bin", type=float, default=0.05, help="Bin width (s) for the wide peri-reward histogram")
    ap.add_argument("--zoom_pre", type=float, default=0.1, help="Seconds before reward in the zoomed raster/histograms")
    ap.add_argument("--zoom_post", type=float, default=0.5, help="Seconds after reward in the zoomed raster/histograms")
    ap.add_argument("--zoom_bin", type=float, default=0.005, help="Bin width (s) for the zoomed histograms")
    ap.add_argument("--test_window", type=float, nargs=2, default=[0.0, 1.0], metavar=("T0", "T1"),
                    help="Post-reward window (s) used for the permutation test and zeroed-fraction")
    ap.add_argument("--lock_pre", type=float, default=0.1,
                    help="Reward-locked window: seconds BEFORE reward onset that still count (default 0.1)")
    ap.add_argument("--lock_post", type=float, default=0.5,
                    help="Reward-locked window: seconds AFTER reward onset that count (default 0.5)")
    ap.add_argument("--n_shuffles", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

    pipeline = load_pipeline_module(args.pipeline_script)
    params = pipeline.parse_params_file(args.params)
    probe_idx = args.probe if args.probe is not None else params["PROBE_INDICES"][0]

    reward_file = args.reward_file or os.path.join(params["OUTPUT_DIR"], "XD1_ON.txt")
    out_dir = args.out_dir or os.path.join(params["OUTPUT_DIR"], "artifact_reward_analysis")
    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, f"imec{probe_idx}_artifacts.csv")
    meta_path = os.path.join(out_dir, f"imec{probe_idx}_artifacts_meta.json")

    # ---- artifacts: detect (slow) or reuse -----------------------------------------------------
    if args.reuse and os.path.isfile(csv_path) and os.path.isfile(meta_path):
        logger.info(f"Reusing cached detection: {csv_path}")
        with open(meta_path) as f:
            meta = json.load(f)
        cached = pd.read_csv(csv_path)
        start_samp = cached["start_sample"].to_numpy(np.int64)
        end_samp = cached["end_sample"].to_numpy(np.int64)
        fs, n_samples = meta["fs"], meta["n_samples"]
    else:
        rec, artifacts = detect_with_pipeline(pipeline, params, probe_idx)
        fs, n_samples = rec.get_sampling_frequency(), rec.get_num_samples(0)
        start_samp, end_samp = artifacts_to_samples(artifacts)
        with open(meta_path, "w") as f:
            json.dump(dict(fs=fs, n_samples=int(n_samples), probe=probe_idx,
                           artifact_proportion=params["ARTIFACT_PROPORTION"],
                           saturation_threshold_uV=params["ARTIFACT_SATURATION_THRESHOLD_UV"]), f, indent=2)

    total_dur = n_samples / fs
    order = np.argsort(start_samp)
    start_samp, end_samp = start_samp[order], end_samp[order]
    start_s, end_s = start_samp / fs, end_samp / fs  # same clock as CatGT event times (0 = first sample)

    # ---- rewards -------------------------------------------------------------------------------
    rewards = np.atleast_1d(np.loadtxt(reward_file, dtype=float)) + args.reward_offset
    rewards = np.sort(rewards)
    logger.info(f"Loaded {len(rewards)} reward times from {reward_file}")
    if len(rewards) and (rewards[-1] > total_dur or rewards[0] < 0):
        logger.warning(
            f"Reward times ({rewards[0]:.1f}-{rewards[-1]:.1f} s) fall outside the recording (0-{total_dur:.1f} s). "
            "Check that they are on the same clock as the .ap file (see --reward_offset)."
        )
        rewards = rewards[(rewards >= 0) & (rewards <= total_dur)]

    # ---- analysis ------------------------------------------------------------------------------
    summary_lines, res = analyse(start_s, end_s, rewards, total_dur, args, probe_idx)
    print("\n".join(summary_lines))

    # per-artifact table
    table = pd.DataFrame(dict(
        start_sample=start_samp, end_sample=end_samp,
        start_s=start_s, end_s=end_s, duration_ms=(end_s - start_s) * 1000,
    ))
    if "t_since_prev" in res:
        table["time_since_prev_reward_s"] = res["t_since_prev"]
        table["time_to_next_reward_s"] = res["t_to_next"]
        table["onset_in_test_window"] = res["in_window"]
        table["signed_latency_to_nearest_reward_s"] = res["near_t"]   # + after reward, - before
        table["reward_locked"] = res["locked"]
    table.to_csv(csv_path, index=False)

    with open(os.path.join(out_dir, f"imec{probe_idx}_artifact_summary.txt"), "w") as f:
        f.write("\n".join(summary_lines) + "\n")

    make_figure(start_s, end_s, rewards, total_dur, res, args, probe_idx, summary_lines,
                os.path.join(out_dir, f"imec{probe_idx}_artifacts_vs_reward.png"),
                os.path.join(out_dir, f"imec{probe_idx}_artifacts_vs_reward.pdf"))
    logger.info(f"Saved results to {out_dir}")


if __name__ == "__main__":
    main()