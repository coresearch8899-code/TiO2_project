"""
Photocatalytic Antibacterial Disinfection -- Figure
=====================================================
2-panel figure: (a) CFU/mL vs. time for S. aureus and E. coli (10^7-dilution
primary series, mean +/- SD of duplicate plates), log-scale y-axis;
(b) bar comparison of the 60-min log10 reduction and pseudo-first-order
inactivation rate constant k_inact between the two organisms.

Uses the same validated categorical palette / style as s05_fig03_05_rsm_surfaces.py and
s01_characterization_figures.py for visual consistency across all
manuscript figures.
"""

import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIGDIR = os.environ.get("FIGDIR", "../figures")
os.makedirs(FIGDIR, exist_ok=True)

BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
MAGENTA, GREEN, VIOLET, RED = "#e87ba4", "#008300", "#4a3aa7", "#e34948"
INK, INK2, MUTED, GRID, BASE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 9,
    "axes.titlesize": 9.5,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.edgecolor": BASE,
    "axes.linewidth": 0.8,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "text.color": INK,
    "axes.labelcolor": INK,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.axisbelow": True,
})

def clean_axes(ax, grid_axis="y"):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, linewidth=0.7, zorder=0, which="both")
    ax.set_axisbelow(True)

def save(fig, name):
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{FIGDIR}/{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {FIGDIR}/{name}.png + .pdf")

with open("../data/antibacterial_results.json") as f:
    R = json.load(f)

ORGS = {"S. aureus": {"color": BLUE, "marker": "o"}, "E. coli": {"color": ORANGE, "marker": "s"}}

fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 3.1), gridspec_kw={"width_ratios": [1.3, 1]})

for org, meta in ORGS.items():
    r = R[org]
    t = np.array(r["time_min"])
    y = np.array(r["cfu_per_mL_mean"])
    yerr = np.array(r["cfu_per_mL_sd"])
    axA.errorbar(t, y, yerr=yerr, color=meta["color"], marker=meta["marker"], ms=4.5,
                 lw=1.3, capsize=2.5, capthick=0.9, elinewidth=0.9, label=org)

axA.set_yscale("log")
axA.set_xlabel("Irradiation time (min)")
axA.set_ylabel("Viable count (CFU mL$^{-1}$)")
axA.set_xlim(-3, 63)
clean_axes(axA, grid_axis="y")
axA.legend(frameon=False, loc="lower left", fontsize=7.6)
axA.text(0.5, 1.04, "(a)", transform=axA.transAxes, fontsize=11, fontweight="bold", ha="center")

orgs = list(ORGS.keys())
x = np.arange(len(orgs))
width = 0.6

log_red = [R[o]["final_log10_reduction_60min"] for o in orgs]
colors = [ORGS[o]["color"] for o in orgs]

axB.bar(x, log_red, width=width, color=colors, zorder=3, edgecolor="white", linewidth=0.6)
for xi, v in zip(x, log_red):
    axB.text(xi, v + 0.006, f"{v:.2f}", ha="center", va="bottom", fontsize=7.8, color=INK)
axB.set_xticks(x)
axB.set_xticklabels([""] * len(orgs))
axB.set_ylabel("60-min log$_{10}$ reduction")
axB.set_ylim(0, max(log_red) * 1.35)
clean_axes(axB, grid_axis="y")
axB.text(0.5, 1.04, "(b)", transform=axB.transAxes, fontsize=11, fontweight="bold", ha="center")

# organism name + k_inact annotation beneath each bar, placed just below the
# x-axis using the axes fraction coordinate system (avoids clipping at y=0)
for xi, o in zip(x, orgs):
    xfrac = (xi - axB.get_xlim()[0]) / (axB.get_xlim()[1] - axB.get_xlim()[0])
    k = R[o]["k_inact_min-1"]
    se = R[o]["k_inact_se"]
    r2 = R[o]["k_inact_R2"]
    axB.text(xfrac, -0.09, o, transform=axB.transAxes, ha="center", va="top",
              fontsize=8, style="italic", color=INK)
    axB.text(xfrac, -0.17, f"k = {k:.4f}±{se:.4f} min$^{{-1}}$\n$R^2$={r2:.2f}",
              transform=axB.transAxes, ha="center", va="top", fontsize=6.4, color=MUTED)

save(fig, "figS1_antibacterial_kinetics")
print("DONE")
