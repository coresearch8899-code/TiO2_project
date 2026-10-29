"""Process-improvement scenarios (revision 2026-10-25), drawn from ../data/lca_scenario_results.json (s13).
(a) GWP per functional unit for each scenario, split into electricity and all other flows;
(b) change relative to the baseline for GWP, the ReCiPe single score, ionizing radiation and terrestrial
    ecotoxicity, the two categories in which rooftop PV electricity raises the burden.
Lettering Arial, at least 8 pt at a print width of 174 mm. Writes ../figures/lca_scenarios.png (+ .pdf).
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

FIGDIR = "../figures"
os.makedirs(FIGDIR, exist_ok=True)
BLUE, ORANGE, AQUA, VIOLET, RED = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7", "#e34948"
INK, INK2, GRID, BASE = "#0b0b0b", "#52514e", "#e1e0d9", "#c3c2b7"
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "axes.edgecolor": BASE, "axes.linewidth": 0.8, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.labelcolor": INK, "figure.facecolor": "white", "savefig.facecolor": "white",
    "axes.axisbelow": True,
})

R = json.load(open("../data/lca_scenario_results.json", encoding="utf-8"))["scenarios"]
labels = [f"{r['id']}  {r['scenario']}" for r in R]
gwp = np.array([r["totals"]["GWP"] for r in R])
gwp_el = np.array([r["electricity_contribution"]["GWP"] for r in R])

fig, (axA, axB) = plt.subplots(1, 2, figsize=(6.85, 3.7), gridspec_kw={"width_ratios": [1.15, 1]})
y = np.arange(len(R))[::-1]
axA.barh(y, gwp_el, height=0.62, color=ORANGE, edgecolor="white", linewidth=0.5, zorder=3, label="Electricity")
axA.barh(y, gwp - gwp_el, left=gwp_el, height=0.62, color=AQUA, edgecolor="white", linewidth=0.5, zorder=3,
         label="All other flows")
for yi, v in zip(y, gwp):
    axA.text(v + 0.12, yi, f"{v:.3f}", va="center", fontsize=8, color=INK)
axA.set_yticks(y)
axA.set_yticklabels(labels, fontsize=8)
axA.set_xlabel("GWP (kg CO$_2$ eq batch$^{-1}$)")
axA.set_xlim(0, gwp.max() * 1.22)
axA.grid(axis="x", color=GRID, linewidth=0.7, zorder=0)
axA.legend(frameon=False, loc="lower left", fontsize=8, bbox_to_anchor=(0.0, 1.0), ncol=2)
axA.set_title("(a)", fontsize=9, fontweight="bold", pad=34)

metrics = [("GWP", "GWP", ORANGE), ("single_score_mPt", "Single score", VIOLET),
           ("IR", "Ionizing radiation", BLUE), ("TET", "Terrestrial ecotoxicity", RED)]
h = 0.19
for k, (key, name, col) in enumerate(metrics):
    vals = np.array([r["pct_change_vs_baseline"][key] for r in R])
    axB.barh(y + (1.5 - k) * h, vals, height=h, color=col, edgecolor="white", linewidth=0.3, zorder=3, label=name)
axB.axvline(0, color=INK2, linewidth=0.8, zorder=4)
axB.set_yticks(y)
axB.set_yticklabels([r["id"] for r in R], fontsize=8)
axB.set_xlabel("Change vs. baseline (%)")
axB.set_xlim(-105, 105)
axB.grid(axis="x", color=GRID, linewidth=0.7, zorder=0)
axB.legend(frameon=False, loc="lower left", fontsize=8, bbox_to_anchor=(0.0, 1.0), ncol=2)
axB.set_title("(b)", fontsize=9, fontweight="bold", pad=34)
for ax in (axA, axB):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{FIGDIR}/lca_scenarios.{ext}", dpi=600, bbox_inches="tight")
print("saved lca_scenarios")
