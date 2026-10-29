"""Graphical abstract for the ET&I submission (revision 2026-10-25): process -> unit operations -> hotspot -> improvement,
with numbers read from ../data/lca_single_process_results.json (s11) and ../data/lca_scenario_results.json (s13).
No atomic Cu-C-TiO2 structure is drawn because copper and carbon incorporation was not determined.
Writes ../figures/graphical_abstract.png (+ .pdf), 2400 x 960 px at 300 dpi.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

FIGDIR = "../figures"
os.makedirs(FIGDIR, exist_ok=True)
ORANGE, AQUA, BLUE, INK, INK2 = "#eb6834", "#1baf7a", "#2a78d6", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 9, "text.color": INK, "figure.facecolor": "white", "savefig.facecolor": "white"})
L = json.load(open("../data/lca_single_process_results.json", encoding="utf-8"))
S = {r["id"]: r for r in json.load(open("../data/lca_scenario_results.json", encoding="utf-8"))["scenarios"]}
sh = L["foreground_inventory"]["energy_shares_pct"]
dry = next(v for k, v in sh.items() if k.startswith("Drying"))
led = next(v for k, v in sh.items() if "LED" in k)

fig, ax = plt.subplots(figsize=(8.0, 3.2))
ax.set_xlim(0, 100); ax.set_ylim(0, 40); ax.axis("off")


def panel(x, title, lines, ec=INK2, fc="white"):
    ax.add_patch(FancyBboxPatch((x, 7.5), 22, 28, boxstyle="round,pad=0.4,rounding_size=1.5", fc=fc, ec=ec, lw=1.1))
    ax.text(x + 11, 32.6, title, ha="center", va="center", fontsize=9.5, fontweight="bold")
    for i, (t, kw) in enumerate(lines):
        ax.text(x + 11, 27.8 - 4.0 * i, t, ha="center", va="center", **kw)


panel(1.0, "Process", [("TiO$_2$ + glucose", {}), ("+ CuSO$_4$", {}), ("↓", {"fontsize": 11}),
                       ("LED/heat treatment", {}), ("60–80 °C, 140 min", {"color": INK2})])
panel(26.0, "Unit operations", [("stirring", {}), ("ultrasound", {}), ("LED treatment", {}),
                                ("centrifugation", {}), ("drying (70 °C, 12 h)", {"color": ORANGE, "fontweight": "bold"})])
panel(51.0, "Hotspot", [("drying", {"fontsize": 9}), (f"{dry:.1f}%", {"fontsize": 12, "color": ORANGE, "fontweight": "bold"}),
                        ("of electricity", {"fontsize": 9, "color": ORANGE}),
                        ("LED treatment", {"fontsize": 9}), (f"{led:.1f}%", {"fontsize": 12, "fontweight": "bold"})],
      ec=ORANGE, fc="#fde3d7")
panel(76.0, "Improvement", [("higher oven loading", {}), ("lower drying energy", {}), ("PV electricity", {}),
                            (f"GWP −{abs(S['S4']['pct_change_vs_baseline']['GWP']):.0f}%", {"color": AQUA, "fontweight": "bold"}),
                            ("(modelled)", {"color": AQUA})],
      ec=AQUA)
for x in (24.2, 49.2, 74.2):
    ax.add_patch(FancyArrowPatch((x, 21.5), (x + 2.0, 21.5), arrowstyle="-|>", mutation_scale=12, lw=1.2, color=INK2))
ax.text(50, 3.0, "The environmental burden was controlled by downstream drying, not by LED treatment.",
        ha="center", fontsize=10, fontweight="bold")
for ext in ("png", "pdf"):
    fig.savefig(f"{FIGDIR}/graphical_abstract.{ext}", dpi=300)
print("saved graphical_abstract")
