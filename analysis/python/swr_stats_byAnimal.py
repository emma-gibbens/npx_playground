"""
Mixed ANOVA pipeline for SWR participation rates
Design: Region (Dorsal/Ventral) as between-subjects
        TrialType (neutral/reward/airpuff) as within-subjects
        Animal as the experimental unit
 
Requires: pandas, numpy, scipy, statsmodels, pingouin
Install:  pip install statsmodels pingouin
"""
import pandas as pd
import numpy as np
from scipy import stats
from itertools import combinations
import statsmodels.formula.api as smf
 
# ── 1. Load and average within animal ────────────────────────────────────────

df = pd.read_csv("D:\\vCA1\\SfNAbstractAnalysis\\SWR_experience_dependence\\SWR_stats_table.csv")

# FIRST A LINEAR MIXED EFFECTS MODEL
# Random intercept for Animal (sessions nested within animals)
model = smf.mixedlm(
    "NormSWR_rate ~ Region * TrialType",
    df,
    groups=df["Session"]
)
result = model.fit()
print(result.summary())

# Now an anova
df_avg = (
    df.groupby(["Animal", "Region", "TrialType"])["NormSWR_rate"]
    .mean()
    .reset_index()
)
 
print("=" * 60)
print("Animal-averaged data")
print("=" * 60)
print(df_avg.to_string(index=False))
print()
 
# Quick sanity check: each animal should appear in only one region
region_per_animal = df_avg.groupby("Animal")["Region"].nunique()
assert (region_per_animal == 1).all(), "Some animals have data in both regions -- check grouping"
 
# ── 2. Descriptive stats ──────────────────────────────────────────────────────
 
print("=" * 60)
print("Descriptive statistics (mean +/- SEM)")
print("=" * 60)
desc = (
    df_avg.groupby(["Region", "TrialType"])["NormSWR_rate"]
    .agg(["mean", "sem", "count"])
    .rename(columns={"mean": "Mean", "sem": "SEM", "count": "N"})
)
print(desc.to_string())
print()
 
# ── 3. Mixed ANOVA ────────────────────────────────────────────────────────────
# Between: Region (Dorsal vs Ventral)
# Within:  TrialType (neutral, reward, airpuff)
# Subject: Animal
 
try:
    import pingouin as pg
 
    print("=" * 60)
    print("Mixed ANOVA (pingouin)")
    print("=" * 60)
    aov = pg.mixed_anova(
        data=df_avg,
        dv="NormSWR_rate",
        within="TrialType",
        between="Region",
        subject="Animal",
        correction=True,        # Greenhouse-Geisser correction for sphericity
    )
    print(aov.to_string(index=False))
    print()
 
    # Sphericity test
    print("=" * 60)
    print("Sphericity test (Mauchly)")
    print("=" * 60)
    spher = pg.sphericity(
        data=df_avg,
        dv="NormSWR_rate",
        within="TrialType",
        subject="Animal",
    )
    print(spher)
    print()
 
except ImportError:
    # Fallback: manual mixed ANOVA via statsmodels AnovaRM for within part
    # and a between-groups comparison
    print("pingouin not found -- running via statsmodels AnovaRM")
    print("(between-subjects Region effect tested separately via t-test)")
    print()
 
    from statsmodels.stats.anova import AnovaRM
 
    print("=" * 60)
    print("Within-subjects ANOVA: TrialType (collapsed across Region)")
    print("=" * 60)
    aov_within = AnovaRM(
        data=df_avg,
        depvar="NormSWR_rate",
        subject="Animal",
        within=["TrialType"],
    ).fit()
    print(aov_within.summary())
    print()
 
    print("=" * 60)
    print("Between-subjects: Region effect (animal means, independent t-test)")
    print("=" * 60)
    animal_means = df_avg.groupby(["Animal", "Region"])["NormSWR_rate"].mean().reset_index()
    dorsal = animal_means[animal_means["Region"] == "Dorsal"]["NormSWR_rate"]
    ventral = animal_means[animal_means["Region"] == "Ventral"]["NormSWR_rate"]
    t, p = stats.ttest_ind(dorsal, ventral)
    print(f"  Dorsal:  mean={dorsal.mean():.4f}, n={len(dorsal)}")
    print(f"  Ventral: mean={ventral.mean():.4f}, n={len(ventral)}")
    print(f"  t = {t:.3f}, p = {p:.4f}")
    print()
 
# ── 4. Post-hoc: TrialType pairwise (within each Region) ─────────────────────
 
print("=" * 60)
print("Post-hoc: TrialType pairwise comparisons (Wilcoxon signed-rank)")
print("Bonferroni correction for 3 comparisons per region")
print("=" * 60)
 
trial_types = ["neutral", "reward", "airpuff"]
pairs = list(combinations(trial_types, 2))
n_comparisons = len(pairs)
 
for region in ["Dorsal", "Ventral"]:
    region_data = df_avg[df_avg["Region"] == region]
    animals = region_data["Animal"].unique()
    print(f"\n  Region: {region} (n={len(animals)} animals)")
    print(f"  {'Comparison':<25} {'W':>6}  {'p (raw)':>8}  {'p (Bonf)':>9}  {'sig':>4}")
    print("  " + "-" * 58)
 
    for tt1, tt2 in pairs:
        vals1 = region_data[region_data["TrialType"] == tt1].set_index("Animal")["NormSWR_rate"]
        vals2 = region_data[region_data["TrialType"] == tt2].set_index("Animal")["NormSWR_rate"]
 
        # Align by animal
        common = vals1.index.intersection(vals2.index)
        v1 = vals1.loc[common].values
        v2 = vals2.loc[common].values
 
        if len(common) < 3:
            print(f"  {tt1} vs {tt2:<15} -- insufficient n ({len(common)})")
            continue
 
        stat, p_raw = stats.wilcoxon(v1, v2, alternative="two-sided")
        p_bonf = min(p_raw * n_comparisons, 1.0)
        sig = "***" if p_bonf < 0.001 else "**" if p_bonf < 0.01 else "*" if p_bonf < 0.05 else "n.s."
        label = f"{tt1} vs {tt2}"
        print(f"  {label:<25} {stat:>6.1f}  {p_raw:>8.4f}  {p_bonf:>9.4f}  {sig:>4}")
 
print()
 
# ── 5. Post-hoc: Region effect within each TrialType ─────────────────────────
 
print("=" * 60)
print("Post-hoc: Region (Dorsal vs Ventral) within each TrialType")
print("Mann-Whitney U (between-subjects, Bonferroni for 3 comparisons)")
print("=" * 60)
 
print(f"\n  {'TrialType':<12} {'U':>6}  {'p (raw)':>8}  {'p (Bonf)':>9}  {'sig':>4}")
print("  " + "-" * 46)
 
for tt in trial_types:
    tt_data = df_avg[df_avg["TrialType"] == tt]
    dorsal = tt_data[tt_data["Region"] == "Dorsal"]["NormSWR_rate"].values
    ventral = tt_data[tt_data["Region"] == "Ventral"]["NormSWR_rate"].values
 
    if len(dorsal) < 2 or len(ventral) < 2:
        print(f"  {tt:<12} -- insufficient n")
        continue
 
    u, p_raw = stats.mannwhitneyu(dorsal, ventral, alternative="two-sided")
    p_bonf = min(p_raw * 3, 1.0)
    sig = "***" if p_bonf < 0.001 else "**" if p_bonf < 0.01 else "*" if p_bonf < 0.05 else "n.s."
    print(f"  {tt:<12} {u:>6.1f}  {p_raw:>8.4f}  {p_bonf:>9.4f}  {sig:>4}")
 
print()
print("Done.")