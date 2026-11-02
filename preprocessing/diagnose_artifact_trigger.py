#!/usr/bin/env python
"""
diagnose_artifact_trigger.py

For every detected artifact, look at the raw whole-probe AP data inside the flagged span (in uV, the
same stream and channels the pipeline's detect_artifacts() uses) and report what actually triggered it:

  * peak |amplitude| (uV) across all channels
  * fraction of channels with |x| >= ARTIFACT_SATURATION_THRESHOLD_UV, at the worst sample
    (the pipeline flags an artifact when this reaches ARTIFACT_PROPORTION, if the rule is amplitude-based)
  * the same, but also requiring the sample to be exactly flat (x[t] == x[t-1]), i.e. railed
  * the same fraction per shank (the detector sees the whole probe, then ONE list is applied to all shanks)
  * how many channels sit at the recording's highest observed value (evidence of a hardware rail)

It also writes the source of SpikeInterface's detect_artifact_periods module next to the results, so you
can read the exact rule your installed version uses.

Usage
-----
    python diagnose_artifact_trigger.py --params D:\\path\\to\\M671_2026_08_27_lfp_ks4_params.txt

Reads the artifact list cached by artifacts_vs_reward.py (imec0_artifacts.csv).
"""
import argparse
import importlib.util
import inspect
import logging
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

PIPELINE_FILENAME = "extract_lfp_and_run_ks4_update.py"
logger = logging.getLogger("diagnose_artifact_trigger")


def load_pipeline_module(path):
    spec = importlib.util.spec_from_file_location("ks4_pipeline", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def open_ap_like_pipeline(params, probe_idx):
    """Whole-probe AP stream, with the same channel selection as the pipeline's detect_artifacts()."""
    import spikeinterface.full as si

    rec = si.read_spikeglx(params["CATGT_OUTPUT_DIR"], stream_name=f"imec{probe_idx}.ap")
    n_total = rec.get_num_channels()
    if n_total >= 1536:
        n_ap = 1536
    elif n_total >= 384:
        n_ap = 384
    else:
        n_ap = n_total
    return rec.select_channels(channel_ids=rec.channel_ids[:n_ap])


def read_rail_uv(catgt_dir, probe_idx, ai_range_override, gain_override):
    """
    Hardware limit of the AP channels in uV, input-referred: imAiRangeMax / imChan0apGain.
    Reads <...>.imec<N>.ap.meta under the CatGT folder, unless both values are given on the command line.
    """
    ai_range = ai_range_override
    gain = gain_override
    if ai_range is None or gain is None:
        meta_path = None
        for folder, _, files in os.walk(catgt_dir):
            for name in files:
                if name.endswith(f"imec{probe_idx}.ap.meta"):
                    meta_path = os.path.join(folder, name)
        if meta_path is None:
            logger.warning("No .ap.meta found; pass --ai_range and --gain to estimate the hardware rail.")
            return None
        with open(meta_path) as f:
            for line in f:
                if "=" not in line:
                    continue
                key, value = line.strip().split("=", 1)
                if key == "imAiRangeMax" and ai_range is None:
                    ai_range = float(value)
                if key == "imChan0apGain" and gain is None:
                    gain = float(value)
        logger.info(f"Read imAiRangeMax={ai_range}, imChan0apGain={gain} from {meta_path}")
    if ai_range is None or gain is None:
        return None
    return ai_range / gain * 1e6


def dump_si_source(out_dir):
    """Save the source of the module that defines detect_artifact_periods so the rule can be read."""
    try:
        import spikeinterface
        import spikeinterface.full as si

        module = inspect.getmodule(si.detect_artifact_periods)
        source = inspect.getsource(module)
        path = os.path.join(out_dir, "si_detect_artifacts_source.py.txt")
        with open(path, "w") as f:
            f.write(f"# spikeinterface {spikeinterface.__version__}\n# {inspect.getsourcefile(module)}\n\n")
            f.write(source)
        logger.info(f"SpikeInterface {spikeinterface.__version__}: wrote detector source to {path}")
        return path
    except Exception as err:
        logger.warning(f"Could not dump SpikeInterface source: {err}")
        return None


def read_scaled(rec, start, end):
    try:
        return rec.get_traces(start_frame=start, end_frame=end, return_scaled=True)
    except Exception:
        logger.warning("Scaled traces unavailable; values are in raw stored units, NOT uV.")
        return rec.get_traces(start_frame=start, end_frame=end)


def analyse(rec_ap, starts, ends, threshold_uv, groups):
    n_samples = rec_ap.get_num_samples()
    n_ch = rec_ap.get_num_channels()
    n_art = len(starts)
    group_values = []
    if groups is not None:
        group_values = sorted(set(groups.tolist()))

    per_channel_max = np.zeros((n_art, n_ch), dtype=np.float32)
    rows = []
    per_group_frac = np.zeros((n_art, len(group_values)))

    for i in range(n_art):
        start = max(1, int(starts[i]))                  # need one earlier sample for the flatness test
        end = min(int(ends[i]), n_samples)
        if end <= start:
            end = start + 1
        traces = read_scaled(rec_ap, start - 1, end)
        x = traces[1:]
        previous = traces[:-1]
        abs_x = np.abs(x)

        over = abs_x >= threshold_uv                    # samples x channels
        flat = x == previous
        n_over_per_sample = over.sum(axis=1)
        n_over_flat_per_sample = (over & flat).sum(axis=1)
        worst = int(np.argmax(n_over_per_sample))

        per_channel_max[i] = abs_x.max(axis=0)

        for g_index, g in enumerate(group_values):
            in_group = groups == g
            per_group_frac[i, g_index] = over[:, in_group].mean(axis=1).max()

        rows.append(dict(
            start_sample=start,
            start_s=None,
            duration_ms=None,
            peak_abs_uv=float(abs_x.max()),
            median_abs_uv_at_worst_sample=float(np.median(abs_x[worst])),
            frac_channels_over_thr_max=float(n_over_per_sample[worst]) / n_ch,
            frac_channels_over_thr_and_flat_max=float(n_over_flat_per_sample.max()) / n_ch,
        ))

    table = pd.DataFrame(rows)

    global_max = float(per_channel_max.max())
    at_rail = np.zeros(n_art, dtype=int)
    for i in range(n_art):
        at_rail[i] = int(np.sum(per_channel_max[i] >= 0.98 * global_max))
    table["n_channels_within_2pct_of_global_peak"] = at_rail
    for g_index, g in enumerate(group_values):
        table[f"frac_over_thr_group{g}"] = per_group_frac[:, g_index]
    return table, global_max, group_values


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description="What triggered each artifact detection?")
    parser.add_argument("--params", required=True)
    parser.add_argument("--pipeline_script", default=os.path.join(here, PIPELINE_FILENAME))
    parser.add_argument("--probe", type=int, default=None)
    parser.add_argument("--artifacts_csv", default=None)
    parser.add_argument("--out_dir", default=None)
    parser.add_argument("--threshold_uv", type=float, default=None, help="Default: ARTIFACT_SATURATION_THRESHOLD_UV")
    parser.add_argument("--ai_range", type=float, default=None, help="imAiRangeMax (V); default: read from .ap.meta")
    parser.add_argument("--gain", type=float, default=None, help="imChan0apGain; default: read from .ap.meta")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

    pipeline = load_pipeline_module(args.pipeline_script)
    params = pipeline.parse_params_file(args.params)
    probe_idx = args.probe if args.probe is not None else params["PROBE_INDICES"][0]
    threshold = args.threshold_uv or params["ARTIFACT_SATURATION_THRESHOLD_UV"]
    proportion = params["ARTIFACT_PROPORTION"]

    out_dir = args.out_dir or os.path.join(params["OUTPUT_DIR"], "artifact_reward_analysis")
    os.makedirs(out_dir, exist_ok=True)
    csv_path = args.artifacts_csv or os.path.join(out_dir, f"imec{probe_idx}_artifacts.csv")
    artifacts = pd.read_csv(csv_path)
    starts = artifacts["start_sample"].to_numpy()
    ends = artifacts["end_sample"].to_numpy()

    rec_ap = open_ap_like_pipeline(params, probe_idx)
    fs = rec_ap.get_sampling_frequency()
    n_ch = rec_ap.get_num_channels()

    groups = rec_ap.get_property("group")
    if groups is None:
        logger.warning("No 'group' (shank) property on the recording; per-shank breakdown skipped.")
    else:
        groups = np.asarray(groups)

    rail_uv = read_rail_uv(params["CATGT_OUTPUT_DIR"], probe_idx, args.ai_range, args.gain)

    dump_si_source(out_dir)

    logger.info(f"Analysing {len(starts)} artifacts on {n_ch} channels ...")
    table, global_max, group_values = analyse(rec_ap, starts, ends, threshold, groups)
    table["start_s"] = table["start_sample"] / fs
    table["duration_ms"] = (ends - starts) * 1000.0 / fs

    try:
        gain = float(rec_ap.get_channel_gains()[0])
    except Exception:
        gain = float("nan")

    n_art = len(table)
    meets = table["frac_channels_over_thr_max"] >= proportion
    lines = []
    lines.append("=" * 78)
    lines.append(f"WHAT TRIGGERED THE DETECTIONS | imec{probe_idx} | {n_art} artifacts | {n_ch} channels")
    lines.append("=" * 78)
    lines.append(f"Pipeline settings: threshold {threshold:g} uV, proportion {proportion:g} "
                 f"(= {proportion * n_ch:.0f} of {n_ch} channels)")
    lines.append(f"Raw gain_to_uV: {gain:.4g} uV per stored unit")
    if rail_uv is not None:
        lines.append(f"Estimated hardware limit (imAiRangeMax / gain): +/- {rail_uv:.0f} uV   "
                     f"-> threshold is {100 * threshold / rail_uv:.0f} % of full scale")
    lines.append("")
    lines.append("Peak |amplitude| inside flagged spans (raw, whole probe, uV)")
    lines.append(f"  median {table['peak_abs_uv'].median():.0f}   "
                 f"min {table['peak_abs_uv'].min():.0f}   max {table['peak_abs_uv'].max():.0f}")
    lines.append(f"  artifacts whose peak never reaches {threshold:g} uV on ANY channel: "
                 f"{int((table['peak_abs_uv'] < threshold).sum())} / {n_art}")
    if rail_uv is not None:
        lines.append(f"  artifacts that reach >= 95% of the hardware limit (true clipping): "
                     f"{int((table['peak_abs_uv'] >= 0.95 * rail_uv).sum())} / {n_art}")
        lines.append(f"  artifacts between the threshold and the limit (large, not clipped): "
                     f"{int(((table['peak_abs_uv'] >= threshold) & (table['peak_abs_uv'] < 0.95 * rail_uv)).sum())} / {n_art}")
    lines.append("")
    lines.append(f"Fraction of channels with |x| >= {threshold:g} uV at the worst sample")
    lines.append(f"  median {table['frac_channels_over_thr_max'].median():.3f}   "
                 f"min {table['frac_channels_over_thr_max'].min():.3f}   "
                 f"max {table['frac_channels_over_thr_max'].max():.3f}")
    lines.append(f"  reach the {proportion:g} proportion: {int(meets.sum())} / {n_art}")
    lines.append(f"  ... and also exactly flat (railed) fraction: median "
                 f"{table['frac_channels_over_thr_and_flat_max'].median():.3f}")
    lines.append("")
    lines.append("Rail evidence")
    lines.append(f"  highest |value| seen on any channel in any flagged span: {global_max:.0f} uV")
    lines.append(f"  channels within 2% of that value, per artifact: median "
                 f"{table['n_channels_within_2pct_of_global_peak'].median():.0f}, "
                 f"max {table['n_channels_within_2pct_of_global_peak'].max()}")
    if len(group_values) > 0:
        lines.append("")
        lines.append("Per shank: median over artifacts of fraction of that shank's channels over threshold")
        only_one = 0
        for i in range(n_art):
            n_groups_over = 0
            for g in group_values:
                if table.loc[i, f"frac_over_thr_group{g}"] >= proportion:
                    n_groups_over += 1
            if n_groups_over == 1:
                only_one += 1
        for g in group_values:
            lines.append(f"  shank {g}: {table[f'frac_over_thr_group{g}'].median():.3f}")
        lines.append(f"  artifacts where only ONE shank has >= {proportion:g} of its channels over threshold: "
                     f"{only_one} / {n_art}")
    lines.append("=" * 78)
    print("\n".join(lines))

    tag = f"imec{probe_idx}"
    table.to_csv(os.path.join(out_dir, f"{tag}_trigger_diagnosis.csv"), index=False)
    with open(os.path.join(out_dir, f"{tag}_trigger_diagnosis.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.2))
    ax = axes[0]
    ax.hist(table["peak_abs_uv"], bins=40, color="0.4")
    ax.axvline(threshold, color="crimson", ls="--", label=f"threshold {threshold:g} uV")
    if rail_uv is not None:
        ax.axvline(rail_uv, color="tab:blue", ls=":", label=f"hardware limit {rail_uv:.0f} uV")
    ax.set_xlabel("Peak |amplitude| in flagged span (uV)")
    ax.set_ylabel("Artifacts")
    ax.legend(frameon=False)
    ax.set_title("A. How large did the raw signal get?", loc="left", fontsize=10)

    ax = axes[1]
    ax.hist(table["frac_channels_over_thr_max"], bins=40, color="0.4")
    ax.axvline(proportion, color="crimson", ls="--", label=f"proportion {proportion:g}")
    ax.set_xlabel(f"Fraction of channels over {threshold:g} uV (worst sample)")
    ax.legend(frameon=False)
    ax.set_title("B. Does the amplitude criterion explain the flags?", loc="left", fontsize=10)

    ax = axes[2]
    if len(group_values) > 0:
        matrix = np.zeros((n_art, len(group_values)))
        for g_index, g in enumerate(group_values):
            matrix[:, g_index] = table[f"frac_over_thr_group{g}"].to_numpy()
        image = ax.imshow(matrix.T, aspect="auto", vmin=0, vmax=1, cmap="magma",
                          extent=[0, n_art, len(group_values) - 0.5, -0.5])
        ax.set_yticks(range(len(group_values)))
        ax.set_yticklabels([f"shank {g}" for g in group_values])
        ax.set_xlabel("Artifact # (in time order)")
        fig.colorbar(image, ax=ax, label="fraction of shank's channels over threshold")
    else:
        ax.axis("off")
    ax.set_title("C. Which shank carries each artifact?", loc="left", fontsize=10)

    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"{tag}_trigger_diagnosis.png"), dpi=170)
    plt.close(fig)
    logger.info(f"Saved results to {out_dir}")


if __name__ == "__main__":
    main()