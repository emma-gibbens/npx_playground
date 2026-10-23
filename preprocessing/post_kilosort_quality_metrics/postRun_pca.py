"""
Finish quality metrics for an existing SpikeInterface analyzer without re-running the pipeline.

Loads the saved analyzer (principal_components and template_metrics already on disk),
computes quality metrics + mahalanobis (isolation_distance, l_ratio), joins template
metrics, then applies the good-unit and noise-flag criteria.

Usage (run in the si_preprocess env):
    python finish_quality_metrics.py "F:\\KS4_outputs\\M671_2026_08_27_imec0_shank0_analyzer"
    python finish_quality_metrics.py <analyzer_folder> --out-dir <folder> --n-jobs 14

Outputs (in --out-dir, default: the analyzer folder):
    quality_metrics_full.csv   all units: quality metrics, template metrics, pass_* columns,
                               good_unit, noise_* columns
    good_unit_ids.csv          just the good unit IDs
"""
import argparse
import os
import pandas as pd
import spikeinterface.full as si

DEFAULT_ANALYZER = r"F:\KS4_outputs\M671_2026_08_27_imec0_shank3_analyzer"

# Keep in sync with the pipeline script.
GOOD_UNIT_CRITERIA = {
    'isi_violations_ratio_max': 0.5,
    'presence_ratio_min': 0.8,
    'snr_min': 1.5,
    'median_amplitude_uv_min': 50.0,
    'amplitude_cutoff_max': 0.1,
}

NOISE_CRITERIA = {
    'peak_to_trough_ms_min': 0.1,
    'peak_to_trough_ms_max': 1.15,
    'num_positive_peaks_max': 2,
    'num_negative_peaks_max': 1,
}

METRIC_NAMES = [
    'num_spikes', 'firing_rate', 'presence_ratio', 'snr',
    'isi_violation', 'rp_violation', 'sliding_rp_violation',
    'amplitude_cutoff', 'amplitude_median', 'amplitude_cv', 'sd_ratio',
]
PC_METRIC_NAMES = ['mahalanobis']   # produces isolation_distance and l_ratio


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('analyzer_folder', nargs='?', default=DEFAULT_ANALYZER)
    parser.add_argument('--out-dir', default=None)
    parser.add_argument('--n-jobs', type=int, default=14)
    args = parser.parse_args()

    print(f"spikeinterface {si.__version__}")
    analyzer = si.load_sorting_analyzer(args.analyzer_folder)
    exts = analyzer.get_loaded_extension_names()
    print(f"Loaded analyzer with {analyzer.get_num_units()} units; extensions: {exts}")

    job_kwargs = dict(n_jobs=args.n_jobs, progress_bar=True)

    # Template metrics: use the saved ones, recompute only if missing
    if 'template_metrics' not in exts:
        print("template_metrics not found, computing...")
        analyzer.compute('template_metrics', verbose=True)

    # Principal components: needed for mahalanobis, recompute only if missing
    if 'principal_components' not in exts:
        print("principal_components not found, computing (slow)...")
        analyzer.compute('principal_components', n_components=5, mode='by_channel_local',
                         verbose=True, **job_kwargs)

    # sd_ratio needs the recording; skip it if the analyzer can't find it on disk
    metric_names = list(METRIC_NAMES)
    if not analyzer.has_recording():
        print("No recording attached to the analyzer; skipping sd_ratio")
        metric_names.remove('sd_ratio')

    analyzer.compute('quality_metrics', metric_names=metric_names + PC_METRIC_NAMES, **job_kwargs)
    qm = analyzer.get_extension('quality_metrics').get_data()
    tm = analyzer.get_extension('template_metrics').get_data()
    qm = qm.join(tm, rsuffix='_tm')

    if 'KSLabel' in analyzer.sorting.get_property_keys():
        qm['ks_label'] = analyzer.sorting.get_property('KSLabel')

    # Good-unit classification
    c = GOOD_UNIT_CRITERIA
    qm['pass_isi'] = qm['isi_violations_ratio'] < c['isi_violations_ratio_max']
    qm['pass_presence'] = qm['presence_ratio'] > c['presence_ratio_min']
    qm['pass_snr'] = qm['snr'] > c['snr_min']
    qm['pass_amplitude'] = qm['amplitude_median'].abs() > c['median_amplitude_uv_min']
    qm['pass_amp_cutoff'] = qm['amplitude_cutoff'] < c['amplitude_cutoff_max']
    pass_cols = ['pass_isi', 'pass_presence', 'pass_snr', 'pass_amplitude', 'pass_amp_cutoff']
    qm['good_unit'] = qm[pass_cols].all(axis=1)

    # Noise flag from waveform shape
    n = NOISE_CRITERIA
    ptv_ms = qm['peak_to_trough_duration'] * 1000.0   # assumes seconds; sanity check below
    med = ptv_ms.median()
    print(f"Median peak_to_trough_duration after *1000: {med:.3f} ms")
    if med > 5:
        print("WARNING: median is implausibly large; peak_to_trough_duration may already be in ms. "
              "Remove the * 1000.0 and rerun.")
    qm['noise_duration'] = (ptv_ms < n['peak_to_trough_ms_min']) | (ptv_ms > n['peak_to_trough_ms_max'])
    qm['noise_n_peaks'] = (
        (qm['num_positive_peaks'] > n['num_positive_peaks_max'])
        | (qm['num_negative_peaks'] > n['num_negative_peaks_max'])
    )
    qm['noise_flag'] = qm['noise_duration'] | qm['noise_n_peaks']

    # Save
    qm.index.name = 'unit_id'
    out_dir = args.out_dir or args.analyzer_folder
    os.makedirs(out_dir, exist_ok=True)
    full_csv = os.path.join(out_dir, 'quality_metrics_full.csv')
    ids_csv = os.path.join(out_dir, 'good_unit_ids.csv')
    qm.to_csv(full_csv)
    qm.index[qm['good_unit']].to_series().to_csv(ids_csv, index=False, header=['unit_id'])

    # Summary
    n_units, n_good = len(qm), int(qm['good_unit'].sum())
    n_noise = int(qm['noise_flag'].sum())
    n_good_noise = int((qm['noise_flag'] & qm['good_unit']).sum())
    print(f"\n{n_good}/{n_units} units pass all good-unit criteria")
    for col in pass_cols:
        print(f"  {col:16s} {int(qm[col].sum()):5d} individually passing")
    print(f"{n_noise}/{n_units} units flagged as noise by waveform shape "
          f"({n_good_noise} also pass good-unit criteria; "
          f"{int(qm['peak_to_trough_duration'].isna().sum())} have NaN shape metrics)")
    print(f"Median isolation distance = {qm['isolation_distance'].median():.1f}, "
          f"median L-ratio = {qm['l_ratio'].median():.3f}")
    print(f"Saved: {full_csv}\nSaved: {ids_csv}")


if __name__ == '__main__':
    main()