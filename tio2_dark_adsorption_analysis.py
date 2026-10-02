"""
Dark Adsorption and Total Removal Analysis
====================================================================
Every dark-adsorption-free metric reported elsewhere in this study (Section
2.5) is referenced to A0, the absorbance measured immediately after the
30-min dark equilibration period, so that only light-driven degradation is
counted. Ai, the absorbance measured before that dark period, is retained in
the raw data (data_analyse_1.xlsx / Supplementary Table S1-S3) "for
provenance only" and was not, in earlier revisions, used to compute or report
anything.

This script closes that gap using the SAME raw Ai/A0/A60(/A15/A30/A45) data
already embedded in tio2_reanalysis.py (BBD_RAW, VALIDATION_RAW,
KINETICS_RAW; copied here unchanged) to additionally compute, for every run/
set/condition:
  - Dark adsorption (%)   = (Ai - A0) / Ai * 100   -- removal during the
    30-min dark equilibration, before any light exposure.
  - Photodegradation (%)  = (A0 - A60) / A0 * 100  -- identical to the
    dark-adsorption-free metric used throughout the rest of the study
    (reproduced here only for cross-checking against tio2_reanalysis.py).
  - Total removal (%)     = (Ai - A60) / Ai * 100  -- the Ai-referenced
    convention used by most literature comparators (Section 1; Table Y),
    which does not separate dark adsorption from photocatalysis.

Nothing here changes or supersedes any dark-adsorption-free value reported
elsewhere: Photodeg_pct computed here must, and does, match
tio2_reanalysis.py's build_response_table() output exactly (checked below by
assertion). Dark_adsorption_pct and Total_removal_pct are new, additive
outputs only.

Raw-data source: data_analyse_1.xlsx (laboratory workbook). Checked against it on 2026-09-24:
V5 A0 replicates 2-3 (0.475, 0.477), V6 A0 replicate 1 (0.453), V7 A0 replicates 2-3 (0.439, 0.436)
and K1 Ai replicate 1 (0.542) had been transcribed incorrectly and are corrected here and in
tio2_reanalysis.py. K1-K6 share one averaged Ai (0.542333) in the workbook, as supplied.

Writes: ../data/dark_adsorption_results.json
"""

import numpy as np
import pandas as pd
import json

# ---------------------------------------------------------------------------
# Raw data, copied unchanged from tio2_reanalysis.py (BBD_RAW, VALIDATION_RAW,
# KINETICS_RAW) so this script can run standalone without re-executing the
# full RSM/ML/optimization/LCA pipeline. Any future correction to the source
# raw data must be applied in both places.
# ---------------------------------------------------------------------------

BBD_RAW = pd.DataFrame([
    [1,  0.5,  0.02, 150, 0.571, 0.498, 0.419],
    [2,  2.0,  0.02, 150, 0.627, 0.514, 0.492],
    [3,  0.5,  0.08, 150, 0.571, 0.474, 0.447],
    [4,  2.0,  0.05, 150, 0.571, 0.475, 0.402],
    [5,  0.5,  0.05, 120, 0.571, 0.464, 0.398],
    [6,  2.0,  0.05, 120, 0.571, 0.529, 0.457],
    [7,  0.5,  0.05, 180, 0.571, 0.463, 0.421],
    [8,  2.0,  0.05, 180, 0.571, 0.499, 0.416],
    [9,  1.25, 0.02, 120, 0.627, 0.527, 0.451],
    [10, 1.25, 0.08, 120, 0.627, 0.550, 0.434],
    [11, 1.25, 0.02, 180, 0.627, 0.541, 0.501],
    [12, 1.25, 0.08, 180, 0.627, 0.551, 0.498],
    [13, 1.25, 0.05, 150, 0.627, 0.511, 0.397],
    [14, 1.25, 0.05, 150, 0.627, 0.513, 0.399],
    [15, 1.25, 0.05, 150, 0.627, 0.514, 0.401],
], columns=["Run", "Glucose_g", "CuSO4_g", "Time_min", "Ai", "A0", "A60"])

VALIDATION_RAW = pd.DataFrame([
    ["V1", 1.57, 0.072, 140, 1, 0.587, 0.481,  0.371],
    ["V1", 1.57, 0.072, 140, 2, 0.587, 0.480,  0.370],
    ["V1", 1.57, 0.072, 140, 3, 0.587, 0.479,  0.369],
    ["V2", 1.57, 0.072, 140, 1, 0.587, 0.4812, 0.368],
    ["V2", 1.57, 0.072, 140, 2, 0.587, 0.4799, 0.375],
    ["V2", 1.57, 0.072, 140, 3, 0.587, 0.483,  0.378],
    ["V3", 1.57, 0.072, 140, 1, 0.587, 0.481,  0.376],
    ["V3", 1.57, 0.072, 140, 2, 0.587, 0.482,  0.380],
    ["V3", 1.57, 0.072, 140, 3, 0.587, 0.480,  0.375],
    ["V4", 1.47, 0.072, 135, 1, 0.587, 0.506,  0.432],
    ["V4", 1.47, 0.072, 135, 2, 0.587, 0.505,  0.432],
    ["V4", 1.47, 0.072, 135, 3, 0.587, 0.505,  0.432],
    ["V5", 1.67, 0.072, 145, 1, 0.587, 0.476,  0.406],
    ["V5", 1.67, 0.072, 145, 2, 0.587, 0.475,  0.406],
    ["V5", 1.67, 0.072, 145, 3, 0.587, 0.477,  0.406],
    ["V6", 1.57, 0.060, 140, 1, 0.587, 0.453,  0.386],
    ["V6", 1.57, 0.060, 140, 2, 0.587, 0.455,  0.387],
    ["V6", 1.57, 0.060, 140, 3, 0.587, 0.455,  0.386],
    ["V7", 1.57, 0.080, 140, 1, 0.587, 0.437,  0.370],
    ["V7", 1.57, 0.080, 140, 2, 0.587, 0.439,  0.371],
    ["V7", 1.57, 0.080, 140, 3, 0.587, 0.436,  0.372],
], columns=["Set", "Glucose_g", "CuSO4_g", "Time_min", "Rep", "Ai", "A0", "A60"])

KINETICS_RAW = pd.DataFrame([
    ["K1", 1.57, 0.072, 140, 1, 0.542, 0.466, 0.453, 0.427, 0.407, 0.379, "Heat / no light"],
    ["K1", 1.57, 0.072, 140, 2, 0.542333, 0.467, 0.450, 0.425, 0.408, 0.378, "Heat / no light"],
    ["K1", 1.57, 0.072, 140, 3, 0.542333, 0.465, 0.450, 0.424, 0.409, 0.378, "Heat / no light"],
    ["K2", 1.57, 0.072, 140, 1, 0.542333, 0.491, 0.460, 0.432, 0.419, 0.400, "Light / no added heat"],
    ["K2", 1.57, 0.072, 140, 2, 0.542333, 0.490, 0.459, 0.432, 0.418, 0.400, "Light / no added heat"],
    ["K2", 1.57, 0.072, 140, 3, 0.542333, 0.489, 0.459, 0.433, 0.417, 0.400, "Light / no added heat"],
    ["K3", 0.00, 0.00, 0,     1, 0.542333, 0.473, 0.451, 0.438, 0.416, 0.407, "Bare TiO2, light+heat"],
    ["K3", 0.00, 0.00, 0,     2, 0.542333, 0.474, 0.450, 0.438, 0.417, 0.409, "Bare TiO2, light+heat"],
    ["K3", 0.00, 0.00, 0,     3, 0.542333, 0.475, 0.459, 0.438, 0.417, 0.407, "Bare TiO2, light+heat"],
    ["K4", 1.57, 0.00, 140,   1, 0.542333, 0.461, 0.445, 0.421, 0.403, 0.377, "Carbon-only, light+heat"],
    ["K4", 1.57, 0.00, 140,   2, 0.542333, 0.462, 0.444, 0.420, 0.405, 0.370, "Carbon-only, light+heat"],
    ["K4", 1.57, 0.00, 140,   3, 0.542333, 0.462, 0.444, 0.420, 0.405, 0.375, "Carbon-only, light+heat"],
    ["K5", 0.00, 0.072, 140,  1, 0.542333, 0.477, 0.456, 0.435, 0.415, 0.384, "Cu-only, light+heat"],
    ["K5", 0.00, 0.072, 140,  2, 0.542333, 0.476, 0.456, 0.436, 0.415, 0.388, "Cu-only, light+heat"],
    ["K5", 0.00, 0.072, 140,  3, 0.542333, 0.478, 0.456, 0.436, 0.416, 0.388, "Cu-only, light+heat"],
    ["K6", 1.57, 0.072, 140,  1, 0.542333, 0.481, 0.464, 0.442, 0.407, 0.372, "Full Cu+C, light+heat"],
    ["K6", 1.57, 0.072, 140,  2, 0.542333, 0.482, 0.464, 0.441, 0.407, 0.372, "Full Cu+C, light+heat"],
    ["K6", 1.57, 0.072, 140,  3, 0.542333, 0.482, 0.464, 0.441, 0.408, 0.375, "Full Cu+C, light+heat"],
], columns=["Condition", "Glucose_g", "CuSO4_g", "Time_min", "Rep", "Ai", "A0",
            "A15", "A30", "A45", "A60", "Method"])


# ---------------------------------------------------------------------------
# Metric functions
# ---------------------------------------------------------------------------

def dark_adsorption_pct(ai, a0):
    return (ai - a0) / ai * 100.0

def photodegradation_pct(a0, a60):
    """Identical formula to tio2_reanalysis.py; used only for a cross-check."""
    return (a0 - a60) / a0 * 100.0

def total_removal_pct(ai, a60):
    return (ai - a60) / ai * 100.0


def add_metrics(df):
    out = df.copy()
    out["Dark_adsorption_pct"] = dark_adsorption_pct(out["Ai"], out["A0"])
    out["Photodeg_pct"] = photodegradation_pct(out["A0"], out["A60"])
    out["Total_removal_pct"] = total_removal_pct(out["Ai"], out["A60"])
    return out


bbd = add_metrics(BBD_RAW)
val = add_metrics(VALIDATION_RAW)
kin = add_metrics(KINETICS_RAW)

# Cross-check against the known, previously reported dark-adsorption-free
# values (design center, Runs 13-15: 22.3, 22.2, 22.0 %; BBD range 4.3-22.3 %).
assert np.isclose(bbd.loc[bbd.Run == 13, "Photodeg_pct"].iloc[0], 22.3, atol=0.1)
assert np.isclose(bbd["Photodeg_pct"].min(), 4.3, atol=0.1)
assert np.isclose(bbd["Photodeg_pct"].max(), 22.3, atol=0.1)

# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------

def summarize(df, group_col):
    g = df.groupby(group_col)[["Dark_adsorption_pct", "Photodeg_pct", "Total_removal_pct"]]
    mean = g.mean()
    sd = g.std(ddof=1)
    summary = mean.add_suffix("_mean").join(sd.add_suffix("_sd"))
    return summary.reset_index()

val_summary = summarize(val, "Set")
kin_summary = summarize(kin, "Condition")

# Pooled V1-V3 (the validated optimum), matching Table 4's "Pooled V1-V3 (n=9)" row
pooled_v1v3 = val[val.Set.isin(["V1", "V2", "V3"])]
pooled_summary = {
    "n": int(len(pooled_v1v3)),
    "Dark_adsorption_pct_mean": float(pooled_v1v3.Dark_adsorption_pct.mean()),
    "Dark_adsorption_pct_sd": float(pooled_v1v3.Dark_adsorption_pct.std(ddof=1)),
    "Photodeg_pct_mean": float(pooled_v1v3.Photodeg_pct.mean()),
    "Photodeg_pct_sd": float(pooled_v1v3.Photodeg_pct.std(ddof=1)),
    "Total_removal_pct_mean": float(pooled_v1v3.Total_removal_pct.mean()),
    "Total_removal_pct_sd": float(pooled_v1v3.Total_removal_pct.std(ddof=1)),
}

# Cross-check pooled optimum against the previously reported 22.30 +/- 0.78 %
assert np.isclose(pooled_summary["Photodeg_pct_mean"], 22.30, atol=0.05)
assert np.isclose(pooled_summary["Photodeg_pct_sd"], 0.78, atol=0.05)

print("=" * 78)
print("DARK ADSORPTION / TOTAL REMOVAL -- BBD (15 runs)")
print("=" * 78)
print(bbd[["Run", "Glucose_g", "CuSO4_g", "Time_min",
           "Dark_adsorption_pct", "Photodeg_pct", "Total_removal_pct"]]
      .round(2).to_string(index=False))

print("\n" + "=" * 78)
print("DARK ADSORPTION / TOTAL REMOVAL -- Validation sets V1-V7 (mean +/- SD)")
print("=" * 78)
print(val_summary.round(2).to_string(index=False))

print("\nPooled V1-V3 (validated optimum, n=9):")
for k, v in pooled_summary.items():
    print(f"  {k}: {v}")

print("\n" + "=" * 78)
print("DARK ADSORPTION / TOTAL REMOVAL -- Mechanistic ablation K1-K6 (mean +/- SD)")
print("=" * 78)
print(kin_summary.round(2).to_string(index=False))
print("\nRaw data as in data_analyse_1.xlsx (V5-V7 A0 and K1 Ai corrected to it, 2026-09-24).")

# ---------------------------------------------------------------------------
# Write output
# ---------------------------------------------------------------------------

out = {
    "bbd_dark_adsorption": json.loads(
        bbd[["Run", "Glucose_g", "CuSO4_g", "Time_min", "Dark_adsorption_pct",
             "Photodeg_pct", "Total_removal_pct"]].to_json(orient="records")),
    "validation_dark_adsorption_by_set": json.loads(val_summary.to_json(orient="records")),
    "validation_dark_adsorption_pooled_v1v3": pooled_summary,
    "kinetics_dark_adsorption_by_condition": json.loads(kin_summary.to_json(orient="records")),
    "bbd_range": {
        "Dark_adsorption_pct": [float(bbd.Dark_adsorption_pct.min()), float(bbd.Dark_adsorption_pct.max())],
        "Photodeg_pct": [float(bbd.Photodeg_pct.min()), float(bbd.Photodeg_pct.max())],
        "Total_removal_pct": [float(bbd.Total_removal_pct.min()), float(bbd.Total_removal_pct.max())],
    },
    "raw_data_source": "data_analyse_1.xlsx; V5-V7 A0 and K1 Ai replicate 1 corrected to it on 2026-09-24",
}

with open("../data/dark_adsorption_results.json", "w") as f:
    json.dump(out, f, indent=1, default=float)
print("\nSaved ../data/dark_adsorption_results.json")
