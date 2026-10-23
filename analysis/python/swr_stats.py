"""
2x3 Mixed ANOVA: Region (between) x TrialType (within)
DV: NormSWR_rate   |   Subject unit: Session
"""

import pandas as pd
import numpy as np
from scipy import stats
import itertools

# df = pd.read_csv("D://vCA1//graveyard//SWR_stats_table.csv")
df = pd.read_csv("D:\\vCA1\\graveyard\\SWR_stats_table_minITI.csv")

# ── Descriptive statistics ─────────────────────────────────────────────────────
print("=" * 65)
print("DESCRIPTIVE STATISTICS")
print("=" * 65)
desc = df.groupby(['Region', 'TrialType'])['NormSWR_rate'].agg(
    N='count', Mean='mean', SD='std', SEM=lambda x: stats.sem(x)
)
print(desc.round(6).to_string())

# ── Build subject-level pivot (subjects × TrialType) per group ─────────────────
# Subject = Session; between factor = Region
pivot = df.pivot_table(index='Session', columns='TrialType',
                       values='NormSWR_rate', aggfunc='mean')
pivot = pivot.join(df[['Session','Region']].drop_duplicates().set_index('Session'))
pivot = pivot.reset_index()

groups = {r: pivot[pivot['Region'] == r] for r in ['Dorsal','Ventral']}
trial_cols = [c for c in pivot.columns if c not in ('Session','Region')]

def mixed_anova(groups, trial_cols):
    """
    Two-way mixed ANOVA (1 between, 1 within).
    Returns a results dict with SS, df, MS, F, p, partial_eta2 for each source.
    """
    all_data = pd.concat(groups.values())
    n_subjects = {k: len(v) for k, v in groups.items()}
    N = sum(n_subjects.values())
    k = len(trial_cols)           # levels of within factor
    a = len(groups)               # levels of between factor
    grand_mean = all_data[trial_cols].values.mean()

    # ── SS_total ──────────────────────────────────────────────────────────────
    SS_total = ((all_data[trial_cols].values - grand_mean)**2).sum()

    # ── SS_between (Region) ───────────────────────────────────────────────────
    SS_A = 0
    for grp_df in groups.values():
        grp_mean = grp_df[trial_cols].values.mean()
        SS_A += grp_df[trial_cols].shape[0] * k * (grp_mean - grand_mean)**2

    # ── SS_subjects_within_groups (error for between factor) ──────────────────
    SS_S_w_A = 0
    for grp_df in groups.values():
        grp_mean = grp_df[trial_cols].values.mean()
        for _, row in grp_df.iterrows():
            subj_mean = row[trial_cols].mean()
            SS_S_w_A += k * (subj_mean - grp_mean)**2

    # ── SS_within (TrialType + interaction + within-error) ────────────────────
    SS_W = SS_total - SS_A - SS_S_w_A  # partitioned further below

    # ── SS_B (TrialType main effect) ──────────────────────────────────────────
    trial_means = all_data[trial_cols].mean()
    SS_B = N * ((trial_means - grand_mean)**2).sum()

    # ── SS_AB (Interaction) ───────────────────────────────────────────────────
    SS_AB = 0
    for grp_df in groups.values():
        grp_mean = grp_df[trial_cols].values.mean()
        cell_means = grp_df[trial_cols].mean()
        for t in trial_cols:
            SS_AB += len(grp_df) * (cell_means[t] - trial_means[t] - grp_mean + grand_mean)**2

    # ── SS_BxS_wA (within-subjects error) ─────────────────────────────────────
    SS_BxS_wA = SS_W - SS_B - SS_AB

    # ── Degrees of freedom ────────────────────────────────────────────────────
    df_A       = a - 1
    df_S_w_A   = N - a
    df_B       = k - 1
    df_AB      = (a - 1) * (k - 1)
    df_BxS_wA  = (k - 1) * (N - a)

    # ── MS and F ──────────────────────────────────────────────────────────────
    MS_A      = SS_A      / df_A
    MS_S_w_A  = SS_S_w_A  / df_S_w_A
    MS_B      = SS_B      / df_B
    MS_AB     = SS_AB     / df_AB
    MS_BxS_wA = SS_BxS_wA / df_BxS_wA

    F_A  = MS_A  / MS_S_w_A
    F_B  = MS_B  / MS_BxS_wA
    F_AB = MS_AB / MS_BxS_wA

    p_A  = stats.f.sf(F_A,  df_A,  df_S_w_A)
    p_B  = stats.f.sf(F_B,  df_B,  df_BxS_wA)
    p_AB = stats.f.sf(F_AB, df_AB, df_BxS_wA)

    # ── Partial eta-squared ───────────────────────────────────────────────────
    np2_A  = SS_A  / (SS_A  + SS_S_w_A)
    np2_B  = SS_B  / (SS_B  + SS_BxS_wA)
    np2_AB = SS_AB / (SS_AB + SS_BxS_wA)

    results = {
        'Region (A)':        dict(SS=SS_A,      df=df_A,      MS=MS_A,      F=F_A,  p=p_A,  np2=np2_A),
        'Error(A)':          dict(SS=SS_S_w_A,  df=df_S_w_A,  MS=MS_S_w_A,  F=None, p=None, np2=None),
        'TrialType (B)':     dict(SS=SS_B,      df=df_B,      MS=MS_B,      F=F_B,  p=p_B,  np2=np2_B),
        'Region×TrialType':  dict(SS=SS_AB,     df=df_AB,     MS=MS_AB,     F=F_AB, p=p_AB, np2=np2_AB),
        'Error(B)':          dict(SS=SS_BxS_wA, df=df_BxS_wA, MS=MS_BxS_wA, F=None, p=None, np2=None),
    }
    return results

results = mixed_anova(groups, trial_cols)

print("\n" + "=" * 65)
print("2×3 MIXED ANOVA")
print("  Between: Region (Dorsal / Ventral)")
print("  Within:  TrialType (neutral / reward / airpuff)")
print("  Subject: Session")
print("=" * 65)
header = f"{'Source':<22} {'SS':>10} {'df':>5} {'MS':>10} {'F':>8} {'p':>8} {'np2':>8}"
print(header)
print("-" * 65)
for src, r in results.items():
    ss = f"{r['SS']:.6f}"
    df_ = str(int(r['df']))
    ms = f"{r['MS']:.6f}"
    F_  = f"{r['F']:.4f}" if r['F'] is not None else "—"
    p_  = f"{r['p']:.4f}" if r['p'] is not None else "—"
    np2 = f"{r['np2']:.4f}" if r['np2'] is not None else "—"
    print(f"{src:<22} {ss:>10} {df_:>5} {ms:>10} {F_:>8} {p_:>8} {np2:>8}")

print("\n  Partial η² conventions: small=0.01, medium=0.06, large=0.14")

# ── Mauchly's Sphericity (Huynh-Feldt approximation) ─────────────────────────
def mauchly_test(groups, trial_cols):
    """Mauchly's W on the pooled within-cell deviation scores."""
    all_data = pd.concat(groups.values())
    X = all_data[trial_cols].values  # N × k matrix
    # Contrast matrix for k-1 orthogonal contrasts
    k = X.shape[1]
    C = np.zeros((k, k-1))
    for j in range(k-1):
        C[j, j]   =  1
        C[j+1:, j] = -1 / (k - j - 1)
    C = C / np.linalg.norm(C, axis=0)
    # Transform to deviation scores
    Y = X @ C   # N × (k-1)
    S = np.cov(Y.T)
    W = np.linalg.det(S) / (np.trace(S)/(k-1))**(k-1)
    N = X.shape[0]
    # Chi-square approx
    df_mauch = k*(k-1)//2 - 1
    chi2 = -(N - 1 - (2*k**2 - 3*k + 3)/(6*(k-1))) * np.log(W)
    p_mauch = stats.chi2.sf(chi2, df_mauch)
    return W, chi2, df_mauch, p_mauch

W, chi2, df_m, p_m = mauchly_test(groups, trial_cols)
print("\n" + "=" * 65)
print("MAUCHLY'S TEST OF SPHERICITY (TrialType)")
print("=" * 65)
print(f"  W = {W:.4f},  χ²({df_m}) = {chi2:.4f},  p = {p_m:.4f}")
if p_m < 0.05:
    print("  ⚠ Sphericity violated — interpret F with Greenhouse-Geisser correction")
else:
    print("  ✓ Sphericity assumption met")

# ── Greenhouse-Geisser epsilon (for reference) ────────────────────────────────
def gg_epsilon(groups, trial_cols):
    all_data = pd.concat(groups.values())
    X = all_data[trial_cols].values
    N, k = X.shape
    row_means = X.mean(axis=1, keepdims=True)
    col_means = X.mean(axis=0, keepdims=True)
    grand = X.mean()
    S = X - row_means - col_means + grand  # residuals
    C = S.T @ S / N
    num = (np.trace(C))**2
    den = (k - 1) * (np.trace(C @ C))
    eps = num / den
    return min(eps, 1.0)

eps_gg = gg_epsilon(groups, trial_cols)
print(f"  Greenhouse-Geisser ε = {eps_gg:.4f}")

# ── POST-HOC TESTS ─────────────────────────────────────────────────────────────

def paired_cohens_d(a, b):
    diff = np.array(a) - np.array(b)
    return diff.mean() / diff.std(ddof=1)

def cohens_d_indep(a, b):
    pooled = np.sqrt(((len(a)-1)*np.std(a,ddof=1)**2 + (len(b)-1)*np.std(b,ddof=1)**2)
                      / (len(a)+len(b)-2))
    return (np.mean(a)-np.mean(b)) / pooled if pooled > 0 else np.nan

# 1. TrialType main effect (paired t-tests, Bonferroni)
print("\n" + "=" * 65)
print("POST-HOC: TrialType main effect")
print("(Paired t-tests across all sessions, Bonferroni corrected for 3 pairs)")
print("=" * 65)
all_pivot = pd.concat(groups.values())
pairs = list(itertools.combinations(trial_cols, 2))
n_pairs = len(pairs)
print(f"{'Comparison':<24} {'t':>7} {'df':>5} {'p-unc':>8} {'p-bonf':>8} {'d':>8}")
print("-" * 65)
for (t1, t2) in pairs:
    a = all_pivot[t1].values
    b = all_pivot[t2].values
    t_s, p_unc = stats.ttest_rel(a, b)
    p_bonf = min(p_unc * n_pairs, 1.0)
    d = paired_cohens_d(a, b)
    n = len(a)
    print(f"{t1:<12} vs {t2:<8} {t_s:>7.4f} {n-1:>5} {p_unc:>8.4f} {p_bonf:>8.4f} {d:>8.4f}")
print("  Cohen's d conventions: small=0.2, medium=0.5, large=0.8")

# 2. Region main effect (independent t-test)
print("\n" + "=" * 65)
print("POST-HOC: Region main effect")
print("(Independent-samples t-test, collapsed across TrialType)")
print("=" * 65)
d_vals = groups['Dorsal'][trial_cols].values.flatten()
v_vals = groups['Ventral'][trial_cols].values.flatten()
t_s, p = stats.ttest_ind(d_vals, v_vals)
d_stat = cohens_d_indep(d_vals, v_vals)
print(f"  Dorsal  mean ± SD: {d_vals.mean():.6f} ± {d_vals.std(ddof=1):.6f}  (N={len(d_vals)})")
print(f"  Ventral mean ± SD: {v_vals.mean():.6f} ± {v_vals.std(ddof=1):.6f}  (N={len(v_vals)})")
print(f"  t({len(d_vals)+len(v_vals)-2}) = {t_s:.4f},  p = {p:.4f},  Cohen's d = {d_stat:.4f}")

# 3. Interaction: TrialType within each Region (paired, Bonferroni)
print("\n" + "=" * 65)
print("POST-HOC: Interaction — TrialType within each Region")
print("(Paired t-tests, Bonferroni corrected for 3 pairs within each region)")
print("=" * 65)
for region, grp_df in groups.items():
    print(f"\n  ── {region} ──")
    print(f"  {'Comparison':<24} {'t':>7} {'df':>5} {'p-unc':>8} {'p-bonf':>8} {'d':>8}")
    print("  " + "-" * 58)
    for (t1, t2) in pairs:
        a = grp_df[t1].values
        b = grp_df[t2].values
        t_s, p_unc = stats.ttest_rel(a, b)
        p_bonf = min(p_unc * n_pairs, 1.0)
        d = paired_cohens_d(a, b)
        n = len(a)
        print(f"  {t1:<12} vs {t2:<8} {t_s:>7.4f} {n-1:>5} {p_unc:>8.4f} {p_bonf:>8.4f} {d:>8.4f}")

# 4. Interaction: Region within each TrialType (independent, Bonferroni)
print("\n" + "=" * 65)
print("POST-HOC: Interaction — Region within each TrialType")
print("(Independent-samples t-tests, Bonferroni corrected for 3 comparisons)")
print("=" * 65)
print(f"{'TrialType':<12} {'t':>7} {'df':>5} {'p-unc':>8} {'p-bonf':>8} {'d':>8}")
print("-" * 55)
for t in trial_cols:
    a = groups['Dorsal'][t].values
    b = groups['Ventral'][t].values
    t_s, p_unc = stats.ttest_ind(a, b)
    p_bonf = min(p_unc * len(trial_cols), 1.0)
    d = cohens_d_indep(a, b)
    df_ = len(a)+len(b)-2
    print(f"{t:<12} {t_s:>7.4f} {df_:>5} {p_unc:>8.4f} {p_bonf:>8.4f} {d:>8.4f}")
    print(f"  {'Dorsal':>8}: {a.mean():.6f} ± {a.std(ddof=1):.6f}  (N={len(a)})")
    print(f"  {'Ventral':>8}: {b.mean():.6f} ± {b.std(ddof=1):.6f}  (N={len(b)})")

print("\n" + "=" * 65)
print("ANALYSIS COMPLETE")
print("=" * 65)
