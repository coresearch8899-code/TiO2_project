"""
Photocatalytic Antibacterial Disinfection Analysis
====================================================
Source data: "New folder" + root uploads -- raw viable-plate-count (CFU) data
for the optimized Cu/carbon-modified TiO2 catalyst (design-center formulation,
1.25 g glucose / 0.05 g CuSO4 / 150 min) against a Gram-positive
(Staphylococcus aureus) and a Gram-negative (Escherichia coli) model bacterium,
under the same 100 W visible-LED irradiation used for the dye-degradation
assay, sampled at t = 0 (dark control), 10, 20, 30, 40, 50, 60 min.

Raw data are serial-dilution spread-plate counts (duplicate plates) at five
ten-fold dilution levels (10^5-10^9). Countable range per standard plate-count
practice is 30-300 CFU/plate; for both organisms the 10^7 dilution falls in
this range (10^5-10^6: too numerous to count, TNTC; 10^8-10^9: too few for
reliable statistics), so the 10^7 dilution is used for all quantitative
analysis below.

CFU/mL = colony count x dilution factor / volume plated. The plated aliquot
volume is not recorded in the source data; a standard 0.1 mL (100 uL)
spread-plate aliquot is assumed throughout (explicitly flagged as an
assumption in the manuscript Methods). This assumption scales every CFU/mL
value by a constant factor and therefore has NO effect on the log10-reduction,
percent-inactivation, or rate-constant results, which depend only on the
count ratios N(t)/N(0).
"""

import numpy as np
import pandas as pd
from scipy import stats
import json

TIMES = np.array([0, 10, 20, 30, 40, 50, 60])
PLATED_VOL_ML = 0.1
DILUTION = 1e7  # countable dilution level used for both organisms

# raw duplicate-plate counts at the 10^7 dilution (primary, within the
# standard 30-300 CFU/plate countable range for both organisms), from the
# source docx. The 10^8 dilution is retained as a secondary corroboration
# series (counts are lower and therefore statistically noisier, but the
# time trend should agree in direction).
RAW = {
    "S. aureus": {
        "rep1": [60, 75, 66, 76, 61, 59, 54],
        "rep2": [92, 73, 68, 57, 71, 66, 61],
    },
    "E. coli": {
        "rep1": [158, 145, 136, 153, 132, 127, 131],
        "rep2": [138, 145, 148, 120, 139, 142, 120],
    },
}
RAW_1E8 = {
    "S. aureus": {"rep1": [9, 6, 5, 7, 5, 4, 4], "rep2": [8, 7, 8, 4, 6, 5, 5]},
    "E. coli": {"rep1": [14, 17, 15, 13, 12, 15, 8], "rep2": [15, 11, 11, 13, 10, 9, 13]},
}

def cfu_per_ml(count):
    return count * DILUTION / PLATED_VOL_ML

results = {}
for organism, reps in RAW.items():
    rep1 = np.array(reps["rep1"], dtype=float)
    rep2 = np.array(reps["rep2"], dtype=float)
    both = np.vstack([rep1, rep2])  # shape (2, 7)
    mean_count = both.mean(axis=0)
    sd_count = both.std(axis=0, ddof=1)

    cfu_mean = cfu_per_ml(mean_count)
    cfu_sd = cfu_per_ml(sd_count)

    N0 = mean_count[0]
    log_reduction = np.log10(N0 / mean_count)  # log10(N0/Nt), per time point
    pct_inactivation = (1 - mean_count / N0) * 100

    # Chick-Watson pseudo-first-order inactivation kinetics:
    # ln(N_t/N0) = -k*t  -->  fit using BOTH replicates as independent y-values
    t_rep = np.concatenate([TIMES, TIMES])
    ln_ratio_rep = np.log(np.concatenate([rep1, rep2]) / N0)
    # exclude t=0 (ln ratio = 0 by definition, uninformative for slope-only fit through origin)
    mask = t_rep > 0
    slope, intercept, r, p, se = stats.linregress(t_rep[mask], ln_ratio_rep[mask])
    k_inact = -slope

    results[organism] = {
        "time_min": TIMES.tolist(),
        "mean_count_10^7dil": mean_count.tolist(),
        "sd_count_10^7dil": sd_count.tolist(),
        "cfu_per_mL_mean": cfu_mean.tolist(),
        "cfu_per_mL_sd": cfu_sd.tolist(),
        "log10_reduction": log_reduction.tolist(),
        "pct_inactivation": pct_inactivation.tolist(),
        "k_inact_min-1": k_inact,
        "k_inact_se": se,
        "k_inact_R2": r ** 2,
        "final_log10_reduction_60min": log_reduction[-1],
        "final_pct_inactivation_60min": pct_inactivation[-1],
        "N0_CFU_per_mL": cfu_mean[0],
    }

print("=" * 78)
print("PHOTOCATALYTIC ANTIBACTERIAL DISINFECTION -- SUMMARY")
print("=" * 78)
for organism, r in results.items():
    print(f"\n--- {organism} ---")
    print(f"  N0 = {r['N0_CFU_per_mL']:.2e} CFU/mL")
    print(f"  60-min log10 reduction = {r['final_log10_reduction_60min']:.2f}")
    print(f"  60-min inactivation = {r['final_pct_inactivation_60min']:.2f}%")
    print(f"  k_inact = {r['k_inact_min-1']:.4f} +/- {r['k_inact_se']:.4f} min^-1, R2 = {r['k_inact_R2']:.4f}")

ratio = results["E. coli"]["k_inact_min-1"] / results["S. aureus"]["k_inact_min-1"]
print(f"\nk_inact(E. coli) / k_inact(S. aureus) = {ratio:.2f}")

# ---------------------------------------------------------------------------
# Secondary corroboration series at the 10^8 dilution: independent evidence
# (lower, noisier counts) used only to check that the direction/magnitude of
# the inactivation trend is not an artifact of the single dilution level
# chosen for the primary analysis. Not used for the headline k_inact values
# reported in the manuscript.
# ---------------------------------------------------------------------------
corroboration = {}
for organism, reps in RAW_1E8.items():
    rep1 = np.array(reps["rep1"], dtype=float)
    rep2 = np.array(reps["rep2"], dtype=float)
    both = np.vstack([rep1, rep2])
    mean_count = both.mean(axis=0)
    N0 = mean_count[0]
    log_reduction = np.log10(N0 / mean_count)
    pct_inactivation = (1 - mean_count / N0) * 100

    t_rep = np.concatenate([TIMES, TIMES])
    ln_ratio_rep = np.log(np.concatenate([rep1, rep2]) / N0)
    mask = t_rep > 0
    slope, intercept, r, p, se = stats.linregress(t_rep[mask], ln_ratio_rep[mask])
    k_inact = -slope

    corroboration[organism] = {
        "mean_count_10^8dil": mean_count.tolist(),
        "log10_reduction": log_reduction.tolist(),
        "pct_inactivation": pct_inactivation.tolist(),
        "k_inact_min-1": k_inact,
        "k_inact_se": se,
        "k_inact_R2": r ** 2,
        "final_log10_reduction_60min": log_reduction[-1],
        "final_pct_inactivation_60min": pct_inactivation[-1],
    }
    results[organism]["corroboration_1e8dil"] = corroboration[organism]

print("\n" + "=" * 78)
print("SECONDARY CORROBORATION SERIES (10^8 dilution) -- DIRECTIONAL CHECK ONLY")
print("=" * 78)
for organism in RAW.keys():
    primary = results[organism]
    corr = corroboration[organism]
    print(f"\n--- {organism} ---")
    print(f"  Primary  (10^7 dil): 60-min log10 reduction = {primary['final_log10_reduction_60min']:.2f}, "
          f"inactivation = {primary['final_pct_inactivation_60min']:.2f}%, "
          f"k_inact = {primary['k_inact_min-1']:.4f} min^-1, R2 = {primary['k_inact_R2']:.4f}")
    print(f"  Corrob.  (10^8 dil): 60-min log10 reduction = {corr['final_log10_reduction_60min']:.2f}, "
          f"inactivation = {corr['final_pct_inactivation_60min']:.2f}%, "
          f"k_inact = {corr['k_inact_min-1']:.4f} min^-1, R2 = {corr['k_inact_R2']:.4f}")
    same_direction = (primary["final_log10_reduction_60min"] > 0) == (corr["final_log10_reduction_60min"] > 0)
    print(f"  Same direction of decline: {same_direction}")

with open("../data/antibacterial_results.json", "w") as f:
    json.dump(results, f, indent=1)
print("\nSaved ../data/antibacterial_results.json")
