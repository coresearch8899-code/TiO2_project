"""Revision figures (2026-09-24), styled as s05_rsm_surface_figures.py, drawn from ../data only.

fig00_graphical_abstract : rebuilt from the current results — full-spectrum LED (380–780 nm; a controlled laboratory
                           surrogate for solar irradiation, not direct sunlight), commercial precursors, validated optimum, single-process LCA; no comparison with
                           a conventional route (none is modeled). All numbers to 3 decimals (author rule).
fig09_mechanistic_ablation: redrawn from kin_summary with TiO₂ and en dashes in the in-figure text.
"""
import json, textwrap
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

M = json.load(open('../data/manuscript_data.json', encoding='utf-8'))
DK = json.load(open('../data/dark_adsorption_results.json', encoding='utf-8'))
L = json.load(open('../data/lca_single_process_results.json', encoding='utf-8'))
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
GREEN, VIOLET, MAGENTA = "#008300", "#4a3aa7", "#e87ba4"
INK, INK2, MUTED, GRID, BASE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "axes.edgecolor": BASE, "axes.linewidth": 0.8, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.labelcolor": INK, "figure.facecolor": "white",
    "savefig.facecolor": "white", "axes.axisbelow": True,
})


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(f"../figures/{name}.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print("saved", name)


# ------------------------------------------------------------------ graphical abstract
pv = DK['validation_dark_adsorption_pooled_v1v3']
pred = M['new_optimum_model_prediction']['Photodeg_pct']
E_TOT = L['foreground_inventory']['Electricity_kWh']
GWP = L['midpoint_total_per_FU']['GWP']
DRY = next(v for k, v in L['foreground_inventory']['energy_shares_pct'].items() if k.startswith('Drying'))

fig, ax = plt.subplots(figsize=(7.4, 6.0))
ax.set_xlim(0, 10); ax.set_ylim(0, 8.1); ax.axis("off")


def flow_box(x, y, w, h, text, fc, fontsize=6.7, fontweight="normal", tc=INK):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.06,rounding_size=0.12",
                                linewidth=1.1, edgecolor=INK, facecolor=fc, zorder=3))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fontsize, color=tc,
            fontweight=fontweight, zorder=4, linespacing=1.45)


ax.text(5.0, 7.95, "Cu/Carbon-Modified TiO₂: Full-Spectrum LED-Assisted Photothermal Synthesis,\n"
        "Constrained Multi-Objective Optimization, and Cradle-to-Gate Life Cycle Assessment",
        fontsize=9.0, color=INK, ha="center", va="top", fontweight="bold", linespacing=1.5)
y_mid, w_box, h_box, gap = 5.15, 2.05, 1.75, 0.28
xs = [0.15 + k * (w_box + gap) for k in range(4)]
flow_box(xs[0], y_mid, w_box, h_box, "Commercial reagents\nTiO₂ + D-glucose\n+ CuSO₄·5H₂O", "#eef2f6")
flow_box(xs[1], y_mid, w_box, h_box, "Full-spectrum\nLED photothermal\nsynthesis, 60–80 °C,\nno calcination\n(Box–Behnken\ndesign)", BLUE, tc="white", fontweight="bold")
flow_box(xs[2], y_mid, w_box, h_box, "Cu/carbon-\nmodified TiO₂\nphotocatalyst", "#eef2f6")
flow_box(xs[3], y_mid, w_box, h_box, "LED photocatalytic\nMB degradation", AQUA, tc="white", fontweight="bold")
for k in range(3):
    ax.add_patch(FancyArrowPatch((xs[k] + w_box + 0.03, y_mid + h_box / 2), (xs[k + 1] - 0.03, y_mid + h_box / 2),
                                 arrowstyle="-|>", mutation_scale=13, linewidth=1.4, color=INK2, zorder=2))
ax.text(0.15, 4.55, "Constrained-optimization condition, experimentally validated\n"
        "(triplicate sets V1–V3, n = 9; inside the convex hull of the Box–Behnken design points):",
        fontsize=7.2, color=INK2, ha="left", va="top", linespacing=1.5)
ax.text(0.15, 3.75, f"Glucose = 1.570 g   |   CuSO₄ = 0.072 g   |   Time = 140.000 min   →   "
        f"Photodegradation = {pv['Photodeg_pct_mean']:.3f} ± {pv['Photodeg_pct_sd']:.3f}%\n"
        f"(60 min, dark adsorption excluded; model-predicted {pred:.3f}%)",
        fontsize=8.1, color=INK, ha="left", va="top", fontweight="bold", linespacing=1.6)
stat_y, stat_w, stat_gap = 1.85, 2.28, 0.18
tiles = [(f"{pv['Photodeg_pct_mean']:.3f}%", "validated photodegradation\n(60 min, n = 9)", ORANGE),
         (f"{E_TOT:.3f} kWh", "process electricity\nper g catalyst", GREEN),
         (f"{GWP:.3f} kg", "CO₂-eq per g catalyst\n(cradle-to-gate GWP)", VIOLET),
         (f"{DRY:.3f}%", "of process electricity\nfrom batch drying", MAGENTA)]
for k, (num, label, col) in enumerate(tiles):
    x0 = 0.15 + k * (stat_w + stat_gap)
    ax.add_patch(FancyBboxPatch((x0, 0.25), stat_w, stat_y, boxstyle="round,pad=0.05,rounding_size=0.1",
                                linewidth=1.0, edgecolor=BASE, facecolor="#fcfcfb", zorder=3))
    ax.text(x0 + stat_w / 2, 0.25 + stat_y * 0.68, num, ha="center", va="center", fontsize=13.5, color=col, fontweight="bold", zorder=4)
    ax.text(x0 + stat_w / 2, 0.25 + stat_y * 0.30, label, ha="center", va="center", fontsize=6.4, color=INK2, zorder=4, linespacing=1.5)
ax.text(5.0, 0.02, "A full-spectrum LED (380–780 nm) was used as a controlled laboratory surrogate for solar irradiation; direct sunlight was not employed.",
        fontsize=6.2, color=INK2, ha="center", va="bottom", style="italic")
save(fig, "fig00_graphical_abstract")

# ------------------------------------------------------------------ Fig. 9 mechanistic ablation
kin = M['kin_summary']
conds = [r['Condition'] for r in kin]
colors = [BLUE if c == "K6" else MUTED for c in conds]
# tick labels describe the treatment actually applied instead of a condition code
ABLATION_TICKS = {"K1": "Heat only (dark)", "K2": "Light only", "K3": "Untreated TiO₂",
                  "K4": "Glucose only", "K5": "CuSO₄ only", "K6": "Glucose + CuSO₄"}
fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 3.55))
xk = np.arange(len(kin))
for ax_, key, lab in ((axA, 'Photodeg_60min_pct', "Photodegradation, 60 min (%)"), (axB, 'k_regression_min-1', "$k_{regression}$ (min$^{-1}$)")):
    ax_.bar(xk, [r[key + '_mean'] for r in kin], yerr=[r[key + '_sd'] for r in kin], color=colors, edgecolor=INK,
            linewidth=0.5, capsize=3, error_kw=dict(elinewidth=0.9, ecolor=INK2), zorder=3)
    ax_.set_xticks(xk); ax_.set_xticklabels([ABLATION_TICKS[c] for c in conds], fontsize=7, rotation=30, ha="right",
                                            rotation_mode="anchor"); ax_.set_ylabel(lab)
    for s in ("top", "right"):
        ax_.spines[s].set_visible(False)
    ax_.grid(axis="y", color=GRID, linewidth=0.7, zorder=0)
axA.set_title("(a)", loc="center", fontsize=9, fontweight="bold")
axB.set_title("(b)", loc="center", fontsize=9, fontweight="bold")
handles = [mpatches.Patch(facecolor=MUTED, edgecolor=INK, label="Control or single-factor treatment"),
           mpatches.Patch(facecolor=BLUE, edgecolor=INK, label="Complete treatment (glucose + CuSO₄)")]
fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=7.3, bbox_to_anchor=(0.5, -0.03))
fig.text(0.5, -0.09, "Treated samples: 1.570 g glucose and/or 0.072 g CuSO₄·5H₂O per 1.000 g TiO₂, 140 min; "
         "untreated TiO₂ did not pass through the treatment step.", ha="center", va="top", fontsize=6.3, color=INK2)
fig.subplots_adjust(bottom=0.34, wspace=0.3)
save(fig, "fig09_mechanistic_ablation")
