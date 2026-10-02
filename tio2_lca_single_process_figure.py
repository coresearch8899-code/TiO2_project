"""
Cradle-to-Gate LCA -- Single-Process Figure (full-spectrum LED-assisted photothermal synthesis)
===============================================================================================
Fig. 10, two panels:
(a) process electricity demand by unit operation (kWh per g catalyst),
    identifying which step actually drives the process's energy footprint;
(b) ReCiPe 2016 endpoint single score decomposed by damage category
    (Human Health / Ecosystems / Resources), identifying which damage
    pathway actually drives the aggregated environmental burden.
Both panels describe the single developed process only; no comparator
route is modeled.
Writes ../figures/fig10_lca_single_process.png (+ .pdf).
"""

import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIGDIR = "../figures"
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
        ax.grid(axis=grid_axis, color=GRID, linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)

def save(fig, name):
    fig.tight_layout()
    fig.savefig(f"{FIGDIR}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{FIGDIR}/{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {FIGDIR}/{name}.png + .pdf")

R = json.load(open("../data/lca_single_process_results.json"))

fig, (axA, axB) = plt.subplots(1, 2, figsize=(8.6, 3.4), gridspec_kw={"width_ratios": [1.25, 1]})

# --- Panel (a): energy by unit operation --------------------------------
energy = R["foreground_inventory"]["energy_breakdown_kWh"]
shares = R["foreground_inventory"]["energy_shares_pct"]
labels_full = list(energy.keys())
short_labels = [
    "Magnetic stirrer", "Ultrasonication", "LED-photothermal\nirradiation",
    "Centrifugation", "Drying oven\n(70 °C, 12 h)",
]
vals = [energy[k] for k in labels_full]
order = np.argsort(vals)
y = np.arange(len(vals))
colors_a = [ORANGE if i == np.argmax(vals) else AQUA for i in range(len(vals))]
axA.barh(y, [vals[i] for i in order], height=0.6, color=[colors_a[i] for i in order],
          zorder=3, edgecolor="white", linewidth=0.6)
axA.set_yticks(y)
axA.set_yticklabels([short_labels[i] for i in order], fontsize=7.6)
axA.set_xlabel("Electricity demand (kWh g$^{-1}$ catalyst)")
for yi, i in zip(y, order):
    axA.text(vals[i] + 0.15, yi, f"{vals[i]:.3f} kWh ({shares[labels_full[i]]:.3f}%)",
              va="center", fontsize=7, color=INK)
axA.set_xlim(0, max(vals) * 1.45)
clean_axes(axA, grid_axis="x")
axA.set_title(r"$\mathbf{(a)}$  Process electricity demand by unit operation",
              fontsize=8.5, color=INK2, loc="left", pad=10)

# --- Panel (b): endpoint single score by damage category ----------------
ss = R["singlescore_total_per_FU_mPt"]
cats = [("mPt_HH", "Human\nHealth", RED), ("mPt_ES", "Ecosystems", GREEN), ("mPt_RS", "Resources", VIOLET)]
totals = [ss[k] for k, _, _ in cats]
total_all = sum(totals)
xb = np.arange(len(cats))
axB.bar(xb, totals, width=0.55, color=[c for _, _, c in cats], zorder=3, edgecolor="white", linewidth=0.6)
axB.set_xticks(xb)
axB.set_xticklabels([lbl for _, lbl, _ in cats], fontsize=8)
axB.set_ylabel("Single score (mPt g$^{-1}$ catalyst)")
for xi, v in zip(xb, totals):
    axB.text(xi, v + total_all * 0.02, f"{v:.3f}\n({100*v/total_all:.3f}%)", ha="center", va="bottom", fontsize=7.4, color=INK)
axB.set_ylim(0, total_all * 1.25)
clean_axes(axB, grid_axis="y")
axB.set_title(r"$\mathbf{(b)}$  Endpoint single score by damage category" + f"\n(total {total_all:.3f} mPt g$^{{-1}}$)",
              fontsize=8.5, color=INK2, loc="left", pad=10)

save(fig, "fig10_lca_single_process")
print("DONE")
