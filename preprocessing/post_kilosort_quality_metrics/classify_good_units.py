"""
Classify "good" units from an existing SpikeInterface quality_metrics/metrics.csv,
using the same criteria as GOOD_UNIT_CRITERIA in extract_lfp_and_run_ks4_update.py.

Usage:
    python classify_good_units.py
    python classify_good_units.py "F:\\KS4_outputs\\<analyzer>\\extensions\\quality_metrics\\metrics.csv"
    python classify_good_units.py <metrics.csv> --out-dir <folder>

Outputs (next to metrics.csv unless --out-dir is given):
    metrics_with_good_units.csv   all units, with pass_* columns and good_unit
    good_unit_ids.csv             just the good unit IDs
"""
import argparse
import os
import pandas as pd

DEFAULT_METRICS_CSV = (
    r"F:\KS4_outputs\M671_2026_08_26_imec0_shank3_analyzer"
    r"\extensions\quality_metrics\metrics.csv"
)

# Keep in sync with GOOD_UNIT_CRITERIA in extract_lfp_and_run_ks4_update.py.
# A unit must pass ALL of these.
GOOD_UNIT_CRITERIA = {
    'isi_violations_ratio_max': 0.5,    # must be < this
    'presence_ratio_min': 0.8,          # must be > this
    'snr_min': 1.5,                     # must be > this
    'median_amplitude_uv_min': 50.0,    # |amplitude_median| must be > this
    'amplitude_cutoff_max': 0.1,        # must be < this
}


def classify_good_units(qm, criteria=GOOD_UNIT_CRITERIA):
    """Add pass_* and good_unit columns. NaN metrics fail their criterion."""
    c = criteria
    qm = qm.copy()
    qm['pass_isi'] = qm['isi_violations_ratio'] < c['isi_violations_ratio_max']
    qm['pass_presence'] = qm['presence_ratio'] > c['presence_ratio_min']
    qm['pass_snr'] = qm['snr'] > c['snr_min']
    qm['pass_amplitude'] = qm['amplitude_median'].abs() > c['median_amplitude_uv_min']
    qm['pass_amp_cutoff'] = qm['amplitude_cutoff'] < c['amplitude_cutoff_max']
    pass_cols = ['pass_isi', 'pass_presence', 'pass_snr', 'pass_amplitude', 'pass_amp_cutoff']
    qm['good_unit'] = qm[pass_cols].all(axis=1)
    return qm, pass_cols


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('metrics_csv', nargs='?', default=DEFAULT_METRICS_CSV)
    parser.add_argument('--out-dir', default=None)
    args = parser.parse_args()

    # First (unnamed) column in SpikeInterface's metrics.csv is the unit id
    qm = pd.read_csv(args.metrics_csv, index_col=0)
    qm.index.name = 'unit_id'

    required = ['isi_violations_ratio', 'presence_ratio', 'snr', 'amplitude_median', 'amplitude_cutoff']
    missing = [m for m in required if m not in qm.columns]
    if missing:
        raise ValueError(f"metrics.csv is missing required columns: {missing}")

    qm, pass_cols = classify_good_units(qm)

    out_dir = args.out_dir or os.path.dirname(os.path.abspath(args.metrics_csv))
    os.makedirs(out_dir, exist_ok=True)
    full_csv = os.path.join(out_dir, 'metrics_with_good_units.csv')
    ids_csv = os.path.join(out_dir, 'good_unit_ids.csv')
    qm.to_csv(full_csv)
    qm.index[qm['good_unit']].to_series().to_csv(ids_csv, index=False, header=['unit_id'])

    n_units, n_good = len(qm), int(qm['good_unit'].sum())
    print(f"{n_good}/{n_units} units pass all good-unit criteria")
    for col in pass_cols:
        print(f"  {col:16s} {int(qm[col].sum()):5d} individually passing")
    print(f"Saved: {full_csv}\nSaved: {ids_csv}")


if __name__ == '__main__':
    main()