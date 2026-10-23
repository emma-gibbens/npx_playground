"""
Read and inspect a SpikeInterface motion-correction output folder.

This matches the folder SpikeInterface writes when you run motion
estimation with a `folder=...` argument (e.g. compute_motion /
correct_motion with preset="medicine"), via Motion.save() plus the
accompanying peaks/parameters/run_times files:

    spikeinterface_info.json
    displacement_seg<N>.npy
    temporal_bins_s_seg<N>.npy
    spatial_bins_um.npy
    parameters.json          (optional)
    run_times.json           (optional)
    peaks.npy                (optional)
    peak_locations.npy       (optional)

The loader below reads these files directly, so it works regardless of
which spikeinterface version wrote them or whether spikeinterface is
even installed in the environment running this script.
"""

import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def load_motion_folder(motion_folder):
    """Load a SpikeInterface Motion output folder into a plain dict."""
    motion_folder = motion_folder+r"\motion"
    motion_folder = Path(motion_folder)

    with open(motion_folder / "spikeinterface_info.json") as f:
        info = json.load(f)
    assert info.get("object") == "Motion", f"{motion_folder} does not look like a Motion folder"

    num_segments = info["num_segments"]
    spatial_bins_um = np.load(motion_folder / "spatial_bins_um.npy")

    displacement, temporal_bins_s = [], []
    for seg in range(num_segments):
        displacement.append(np.load(motion_folder / f"displacement_seg{seg}.npy"))
        temporal_bins_s.append(np.load(motion_folder / f"temporal_bins_s_seg{seg}.npy"))

    result = {
        "info": info,
        "direction": info["direction"],
        "interpolation_method": info["interpolation_method"],
        "num_segments": num_segments,
        "spatial_bins_um": spatial_bins_um,
        "temporal_bins_s": temporal_bins_s,
        "displacement": displacement,
        "rigid": spatial_bins_um.shape[0] == 1,
    }

    # Optional extras saved alongside the Motion object itself
    params_path = motion_folder / "parameters.json"
    if params_path.exists():
        with open(params_path) as f:
            result["parameters"] = json.load(f)

    run_times_path = motion_folder / "run_times.json"
    if run_times_path.exists():
        with open(run_times_path) as f:
            result["run_times"] = json.load(f)

    peaks_path = motion_folder / "peaks.npy"
    if peaks_path.exists():
        result["peaks"] = np.load(peaks_path)

    peak_locations_path = motion_folder / "peak_locations.npy"
    if peak_locations_path.exists():
        result["peak_locations"] = np.load(peak_locations_path)

    return result


def try_load_as_spikeinterface_motion(motion_folder):
    """Optional: wrap the folder in SpikeInterface's own Motion object/
    motion_info dict, which gives you its built-in interpolation and
    plotting helpers. Returns None if spikeinterface isn't installed or
    the loader API in your version differs.
    """
    try:
        import spikeinterface.full as si
        return si.load_motion_info(motion_folder)
    except Exception:
        pass

    try:
        from spikeinterface.sortingcomponents.motion import Motion
        return Motion.load(motion_folder)
    except Exception:
        return None


def plot_motion(motion_data, title=None):
    disp = motion_data["displacement"][0]   # (n_time_bins, n_spatial_bins)
    t = motion_data["temporal_bins_s"][0]    # (n_time_bins,)
    spatial_bins = motion_data["spatial_bins_um"]

    fig, ax = plt.subplots(figsize=(9, 4))

    if motion_data["rigid"]:
        ax.plot(t, disp[:, 0], lw=1)
        ax.set_ylabel(f"displacement ({motion_data['direction']}, \u00b5m)")
    else:
        im = ax.imshow(
            disp.T,
            aspect="auto",
            origin="lower",
            extent=[t[0], t[-1], spatial_bins[0], spatial_bins[-1]],
            cmap="RdBu_r",
        )
        ax.set_ylabel(f"depth ({motion_data['direction']}, \u00b5m)")
        fig.colorbar(im, ax=ax, label="displacement (\u00b5m)")

    ax.set_xlabel("time (s)")
    ax.set_title(title or "Estimated motion")
    fig.tight_layout()
    return fig, ax


def build_motion_summary(motion_data, source_folder_name=None):
    """Build a JSON-serializable summary dict of the motion estimate
    (metadata, per-segment displacement stats, estimation parameters,
    run times, and peak count).
    """
    summary = {
        "source_folder": source_folder_name,
        "direction": motion_data["direction"],
        "interpolation_method": motion_data["interpolation_method"],
        "num_segments": motion_data["num_segments"],
        "rigid": motion_data["rigid"],
        "n_spatial_bins": int(motion_data["spatial_bins_um"].shape[0]),
        "segments": [],
    }

    for seg, (t, d) in enumerate(zip(motion_data["temporal_bins_s"], motion_data["displacement"])):
        summary["segments"].append(
            {
                "segment_index": seg,
                "n_time_bins": int(d.shape[0]),
                "start_s": float(t[0]),
                "end_s": float(t[-1]),
                "displacement_min_um": float(d.min()),
                "displacement_max_um": float(d.max()),
                "displacement_mean_um": float(d.mean()),
                "displacement_std_um": float(d.std()),
            }
        )

    if "parameters" in motion_data:
        params = motion_data["parameters"]
        summary["estimate_motion_kwargs"] = params.get("estimate_motion_kwargs", {})
        summary["detect_kwargs"] = params.get("detect_kwargs", {})
        summary["localize_peaks_kwargs"] = params.get("localize_peaks_kwargs", {})
        summary["sampling_frequency"] = params.get("sampling_frequency")

    if "run_times" in motion_data:
        summary["run_times_s"] = motion_data["run_times"]

    if "peaks" in motion_data:
        summary["n_peaks"] = int(motion_data["peaks"].shape[0])

    return summary


def save_motion_outputs(motion_data, motion_folder, output_dir=None, fig=None):
    """Save everything needed to revisit this motion estimate without
    re-running anything: a JSON summary, the displacement trace as a
    plain CSV (opens in Excel or any text editor), and the plot as a
    PNG. Defaults to a `motion_summary` subfolder right next to the
    source data so it's easy to find later.
    """
    motion_folder = Path(motion_folder)
    output_dir = Path(output_dir) if output_dir else motion_folder / "motion_summary"
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) Summary — readable as plain text, parseable as JSON
    summary = build_motion_summary(motion_data, source_folder_name=motion_folder.name)
    summary_path = output_dir / "motion_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    # 2) Displacement trace as CSV
    t = motion_data["temporal_bins_s"][0]
    d = motion_data["displacement"][0]
    if motion_data["rigid"]:
        header = "time_s,displacement_um"
        data = np.column_stack([t, d[:, 0]])
    else:
        spatial_bins = motion_data["spatial_bins_um"]
        header = "time_s," + ",".join(f"depth_{int(b)}um" for b in spatial_bins)
        data = np.column_stack([t, d])
    csv_path = output_dir / "displacement_trace.csv"
    np.savetxt(csv_path, data, delimiter=",", header=header, comments="", fmt="%.6f")

    # 3) The plot
    png_path = output_dir / "motion_plot.png"
    if fig is not None:
        fig.savefig(png_path, dpi=150, bbox_inches="tight")

    print(f"\nSaved summary, displacement CSV, and plot to:\n  {output_dir}")
    return {"summary": summary_path, "csv": csv_path, "plot": png_path}


if __name__ == "__main__":
    motion_folder = r"F:\KS4_outputs\medicine\M692_2026_08_17_imec0_shank3_motion_corrected"

    motion_data = load_motion_folder(motion_folder)

    print(f"direction:              {motion_data['direction']}")
    print(f"interpolation method:   {motion_data['interpolation_method']}")
    print(f"segments:               {motion_data['num_segments']}")
    print(f"rigid:                  {motion_data['rigid']}")
    print(f"spatial bins:           {motion_data['spatial_bins_um'].shape[0]}")
    for seg, (t, d) in enumerate(zip(motion_data["temporal_bins_s"], motion_data["displacement"])):
        print(
            f"  segment {seg}: {d.shape[0]} time bins, "
            f"{t[0]:.1f}-{t[-1]:.1f} s, "
            f"displacement range {d.min():.2f} to {d.max():.2f} \u00b5m"
        )

    if "parameters" in motion_data:
        est = motion_data["parameters"].get("estimate_motion_kwargs", {})
        print(f"estimation method:      {est.get('method')}, rigid={est.get('rigid')}")

    if "run_times" in motion_data:
        for step, seconds in motion_data["run_times"].items():
            print(f"  {step}: {seconds:.1f} s")

    if "peaks" in motion_data:
        print(f"peaks used:             {motion_data['peaks'].shape[0]:,}")

    fig, ax = plot_motion(motion_data, title=Path(motion_folder).name)
    save_motion_outputs(motion_data, motion_folder, fig=fig)
    plt.show()

    # Optional: also load it as SpikeInterface's native Motion object
    # (useful if you want to call its own interpolation/plotting methods)
    # si_motion_info = try_load_as_spikeinterface_motion(motion_folder)