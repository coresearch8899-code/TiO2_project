"""Cradle-to-gate system boundary and process flow (revision 2026-10-25), with the foreground inventory per
functional unit (1 g catalyst) read from ../data/lca_single_process_results.json (s11).
Lettering Arial, at least 8 pt at a print width of 174 mm. Writes ../figures/lca_system_boundary.png (+ .pdf).
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

FIGDIR = "../figures"
os.makedirs(FIGDIR, exist_ok=True)
ORANGE, AQUA, BLUE, INK, INK2, BASE = "#eb6834", "#1baf7a", "#2a78d6", "#0b0b0b", "#52514e", "#c3c2b7"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 8, "text.color": INK, "figure.facecolor": "white", "savefig.facecolor": "white"})

R = json.load(open("../data/lca_single_process_results.json", encoding="utf-8"))
inv = R["foreground_inventory"]
E = list(inv["energy_breakdown_kWh"].values())  # stirring, ultrasonication, LED, centrifugation, drying

fig, ax = plt.subplots(figsize=(6.85, 3.3))
ax.set_xlim(0, 100); ax.set_ylim(0, 44); ax.axis("off")


def box(x, y, w, h, text, fc="white", ec=INK2, lw=0.9, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2", fc=fc, ec=ec, lw=lw))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=8, fontweight="bold" if bold else "normal")


def arrow(x0, y0, x1, y1, color=INK2):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=8, lw=0.9, color=color))


ax.add_patch(Rectangle((1, 6), 98, 37, fill=False, ls=(0, (4, 3)), lw=1.0, ec=INK2))
ax.text(2.5, 40.6, "Cradle-to-gate system boundary (FU = one batch, 1.000 g TiO$_2$ input)", fontsize=8, fontweight="bold", color=INK2)

steps = [("Suspension\npreparation\n(stirring, 30 min)", E[0]), ("Ultrasonication\n(30 min)", E[1]),
         ("LED photothermal\ntreatment (140 min,\n60–80 °C)", E[2]), ("Centrifugation\nand washing\n(10 min)", E[3]),
         ("Drying\n(70 °C, 12 h)", E[4])]
x0, w, gap, y = 3.5, 16.2, 3.0, 17
for i, (name, kwh) in enumerate(steps):
    x = x0 + i * (w + gap)
    hot = i == 4
    box(x, y, w, 11, name, fc="#fde3d7" if hot else "white", ec=ORANGE if hot else INK2, lw=1.2 if hot else 0.9)
    ax.text(x + w / 2, 32.8, f"{kwh:.3f} kWh", ha="center", fontsize=8, color=ORANGE if hot else INK2,
            fontweight="bold" if hot else "normal")
    arrow(x + w / 2, 31.8, x + w / 2, 28.6, ORANGE if hot else INK2)
    if i < 4:
        arrow(x + w + 0.5, y + 5.5, x + w + gap - 0.5, y + 5.5)
ax.text(50, 36.6, f"Thailand grid electricity, low voltage; rated-power upper bound {inv['Electricity_kWh']:.3f} kWh per batch", ha="center",
        fontsize=8, color=INK2)

ax.text(4, 9.2, f"Inputs: TiO$_2$ {inv['TiO2_g']:.3f} g; D-glucose {inv['Glucose_g']:.3f} g; "
        f"CuSO$_4$·5H$_2$O {inv['CuSO4_5H2O_g']:.3f} g; deionized water {inv['DI_water_kg']:.3f} L; "
        "transport and packaging", fontsize=8)
ax.text(4, 12.4, f"Outputs: precursor-treated TiO$_2$ (gate; recovered mass not measured); wastewater "
        f"{inv['Wastewater_m3'] * 1000:.3f} L to treatment (inside boundary)", fontsize=8)
ax.text(50, 2.0, "Outside the boundary: equipment manufacture, photocatalytic use, catalyst reuse and end-of-life",
        ha="center", fontsize=8, color=INK2, style="italic")
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{FIGDIR}/lca_system_boundary.{ext}", dpi=600, bbox_inches="tight")
print("saved lca_system_boundary")
