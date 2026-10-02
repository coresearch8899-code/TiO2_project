"""
Mechanistic Ablation -- Formal Additive-Model (Bliss-Type) Interaction Test
=============================================================================
Tests whether the complete Cu-and-carbon system (K6) exceeds the bare-TiO2
baseline (K3) by MORE than the sum of the two single-factor excesses
(K4-K3 for carbon alone, K5-K3 for Cu alone), which is the standard
operational definition of a super-additive ("synergistic") interaction
between two factors added to a common baseline.

Interaction term:  I = (K6 - K3) - [(K4 - K3) + (K5 - K3)] = K6 - K4 - K5 + K3

I > 0, statistically significant  -> super-additive (synergistic) interaction
I not significantly different from 0 -> consistent with a purely additive
                                          (independent, non-interacting) combination
I < 0, statistically significant  -> sub-additive (antagonistic) interaction

This directly tests the claim -- made in an earlier draft of the manuscript --
that the mechanistic ablation "confirmed a genuine synergy" between copper and
carbon. A rigorous test is required because K6 exceeding every single-factor
control (K4, K5) individually is necessary but NOT sufficient evidence of a
super-additive interaction; it is equally consistent with two independent,
purely additive contributions.

Variance of I is estimated from the four condition-level triplicate SDs
(K3, K4, K5, K6; n = 3 each) assuming independence between conditions, and
significance is assessed with a Welch-Satterthwaite-corrected t-test (unequal
variances, unequal effective df) since the four triplicate SDs differ
considerably (photodegradation SD ranges from 0.30% for K3 to 0.86% for K4).
"""

import json
import numpy as np
from scipy import stats

D = json.load(open("../data/manuscript_data.json", encoding="utf-8"))
kin = {r["Condition"]: r for r in D["kin_summary"]}
K3, K4, K5, K6 = kin["K3"], kin["K4"], kin["K5"], kin["K6"]
N = 3


def welch_satterthwaite_df(sd_n_pairs):
    num = sum(sd**2 / n for sd, n in sd_n_pairs) ** 2
    den = sum((sd**2 / n) ** 2 / (n - 1) for sd, n in sd_n_pairs)
    return num / den


def interaction_test(mean_key, sd_key, label, unit=""):
    interaction = K6[mean_key] - K4[mean_key] - K5[mean_key] + K3[mean_key]
    var_sum = sum((r[sd_key] ** 2) / N for r in (K6, K4, K5, K3))
    se = np.sqrt(var_sum)
    t = interaction / se
    df = welch_satterthwaite_df([(r[sd_key], N) for r in (K6, K4, K5, K3)])
    p = 2 * (1 - stats.t.cdf(abs(t), df))
    carbon_excess = K4[mean_key] - K3[mean_key]
    cu_excess = K5[mean_key] - K3[mean_key]
    combined_excess = K6[mean_key] - K3[mean_key]
    print(f"--- {label} ---")
    print(f"  K3 (baseline)     = {K3[mean_key]:.5f} {unit}")
    print(f"  K4 (carbon-only), excess over K3 = {carbon_excess:.5f} {unit}")
    print(f"  K5 (Cu-only),     excess over K3 = {cu_excess:.5f} {unit}")
    print(f"  K6 (full),        excess over K3 = {combined_excess:.5f} {unit}  (observed combined effect)")
    print(f"  Sum of individual excesses (additive-model expectation) = {carbon_excess + cu_excess:.5f} {unit}")
    print(f"  Interaction term I = K6 - K4 - K5 + K3 = {interaction:.5f} {unit}")
    print(f"  SE(I) = {se:.5f}, t = {t:.3f}, Welch-Satterthwaite df = {df:.2f}, two-sided p = {p:.4f}")
    verdict = ("statistically significant super-additive (synergistic) interaction" if (p < 0.05 and interaction > 0)
               else "statistically significant sub-additive (antagonistic) interaction" if (p < 0.05 and interaction < 0)
               else "not statistically distinguishable from a purely additive (non-interacting) combination")
    print(f"  Verdict: {verdict}\n")
    return {
        "baseline_K3": K3[mean_key], "carbon_excess": carbon_excess, "cu_excess": cu_excess,
        "combined_excess": combined_excess, "additive_expectation": carbon_excess + cu_excess,
        "interaction_term": interaction, "se": se, "t": t, "df": df, "p_value": p, "verdict": verdict,
    }


print("=" * 78)
print("MECHANISTIC ABLATION: FORMAL ADDITIVE-MODEL INTERACTION TEST (K1-K6)")
print("=" * 78)
results = {
    "photodegradation_pct": interaction_test("Photodeg_60min_pct_mean", "Photodeg_60min_pct_sd", "Photodegradation efficiency, 60 min", "%"),
    "k_regression_min-1": interaction_test("k_regression_min-1_mean", "k_regression_min-1_sd", "Pseudo-first-order rate constant", "min^-1"),
}

with open("../data/synergy_analysis_results.json", "w") as f:
    json.dump(results, f, indent=1)
print("Saved ../data/synergy_analysis_results.json")
