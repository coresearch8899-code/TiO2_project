"""Fig. 8 (formerly Fig. 6): validation sets V1–V7, redrawn from manuscript_data.json with the same styling as
s05_rsm_surface_figures.py (FIGURE 5 block), fixing the in-figure text (fix-list m29): CuSO₄ subscript, true minus signs,
en dash, and 3-decimal reference value (author instruction 2026-09-24)."""
import json, textwrap
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

M = json.load(open('../data/manuscript_data.json', encoding='utf-8'))
BLUE, ORANGE, VIOLET = "#2a78d6", "#eb6834", "#4a3aa7"
INK, INK2, GRID, BASE = "#0b0b0b", "#52514e", "#e1e0d9", "#c3c2b7"
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "axes.edgecolor": BASE, "axes.linewidth": 0.8, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.labelcolor": INK, "figure.facecolor": "white",
    "savefig.facecolor": "white", "axes.axisbelow": True, "axes.unicode_minus": True,
})

val = M['val_summary']
center = M['center_photodeg_mean']
sets = [r['Set'] for r in val]
fig, ax = plt.subplots(figsize=(7.2, 4.1))
xv = np.arange(len(val))
ax.bar(xv, [r['Photodeg_mean'] for r in val], yerr=[r['Photodeg_sd'] for r in val],
       color=[BLUE if s in ("V1", "V2", "V3") else ORANGE for s in sets], edgecolor=INK, linewidth=0.5,
       capsize=3, error_kw=dict(elinewidth=0.9, ecolor=INK2), zorder=3, width=0.55)
ax.scatter(xv, [r['Model_prediction_at_point_pct'] for r in val], marker="D", s=42, color=VIOLET, edgecolor=INK,
           linewidth=0.7, zorder=4, label="RSM prediction at that set’s own condition")
ax.axhline(center, color=INK, linestyle="--", linewidth=1.1, zorder=2,
           label=f"Original design center ({center:.3f}%, n = 3)")
# tick labels give the tested condition itself (glucose, CuSO4, time), not a set code
REPLICATE = {"V1": 1, "V2": 2, "V3": 3}
ticks = [f"{r['Glucose_g']:.3f} g\n{r['CuSO4_g']:.3f} g\n{r['Time_min']:.0f} min"
         + (f"\nreplicate {REPLICATE[r['Set']]}" if r['Set'] in REPLICATE else "") for r in val]
ax.set_xticks(xv); ax.set_xticklabels(ticks, fontsize=7.2)
ax.set_ylim(0, 29)
ax.set_ylabel("Photodegradation, 60 min (%)")
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.grid(axis="y", color=GRID, linewidth=0.7, zorder=0)
ax.legend(frameon=False, fontsize=6.9, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.0))
fig.text(0.5, 0.005, "Tick labels: glucose, CuSO₄·5H₂O and LED-treatment time per 1.000 g TiO₂. "
         "Blue: selected condition; orange: perturbed conditions.", ha="center", va="bottom", fontsize=6.6, color=INK2)
fig.subplots_adjust(bottom=0.27)
for ext in ("png", "pdf"):
    fig.savefig(f"../figures/fig08_validation.{ext}", dpi=600, bbox_inches="tight")
print("saved fig08_validation")
