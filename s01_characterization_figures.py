"""
TiO2 Cu/Carbon-Modified Photocatalyst -- Materials Characterization Figures
============================================================================
Builds the two new characterization figures added to the manuscript:

  fig0b_characterization_xrd_ftir_sem.png  -- 4-panel: (a) XRD, (b) FTIR,
                                               (c) SEM 15kx, (d) SEM 30kx
  fig0c_drs_tauc.png                        -- 2-panel: (a) UV-Vis DRS
                                               reflectance, (b) Tauc plot
                                               (indirect transition) with
                                               band-gap extrapolation

Source data: New folder/ (staged from the user's Windows machine)
  - 496398.png              : XRD patterns, runs N7-N12, vs JCPDS 00-004-0477
  - 496397.png              : FTIR spectra, runs N7-N12
  - SEM TiO2(1)/.../SEM/*.bmp : SEM micrographs, run 13 (15 kx & 30 kx)
  - TiO2(2)/TiO2/%R.xlsx     : raw %R spectra + Kubelka-Munk / Tauc columns
    per BBD run (sheets 'TiO2-1' ... 'TiO2-15')

Uses the same validated categorical palette / style as s05_fig03_05_rsm_surfaces.py so
all manuscript figures are visually consistent.
"""

import os
import numpy as np
import pandas as pd
import openpyxl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from PIL import Image

FIGDIR = os.environ.get("FIGDIR", "../figures")
os.makedirs(FIGDIR, exist_ok=True)
SRC = "../data/"   # revision package layout (author-supplied files, 2026-09-24)

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

def save(fig, name, tight=True):
    if tight:
        fig.tight_layout()
    fig.savefig(f"{FIGDIR}/{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{FIGDIR}/{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {FIGDIR}/{name}.png + .pdf")

print("=" * 78)
print("GENERATING MATERIALS-CHARACTERIZATION FIGURES")
print("=" * 78)

# ---------------------------------------------------------------------------
# Panel-labelled composite: XRD (source photo) + FTIR (source photo)
#                           + SEM 15kx + SEM 30kx
# ---------------------------------------------------------------------------
xrd_img = mpimg.imread(SRC + "496398_0.png")
ftir_img = mpimg.imread(SRC + "496397_0.png")
sem15_img = np.asarray(Image.open(
    SRC + "SEM TiO2/69.06.09 (IT800-0128)/SEM/TiO2-No.13-15kx-1.bmp").convert("RGB"))
sem30_img = np.asarray(Image.open(
    SRC + "SEM TiO2/69.06.09 (IT800-0128)/SEM/TiO2-No.13-30kx-1.bmp").convert("RGB"))

fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.6))
panel_imgs = [xrd_img, ftir_img, sem15_img, sem30_img]
panel_labels = ["a", "b", "c", "d"]
for ax, img, lab in zip(axes.flat, panel_imgs, panel_labels):
    ax.imshow(img)
    ax.axis("off")
    ax.text(0.5, 0.98, f"({lab})", transform=ax.transAxes, fontsize=11, fontweight="bold",
            va="top", ha="center", color=INK,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75))
fig.subplots_adjust(wspace=0.03, hspace=0.03)
save(fig, "fig01_characterization_xrd_ftir_sem", tight=False)

# ---------------------------------------------------------------------------
# UV-Vis DRS reflectance + Tauc (indirect transition) band-gap extrapolation
# for three representative Box-Behnken runs
# ---------------------------------------------------------------------------
REPS = {
    2:  {"label": "Run 2 (low activity)",  "color": ORANGE},
    10: {"label": "Run 10 (high CuSO$_4$)", "color": AQUA},
    13: {"label": "Run 13 (design centre)", "color": BLUE},
}
# Table values (DRS onset method) for consistency with the reported Eg table
EG_INDIRECT = {2: 3.23, 10: 3.21, 13: 3.22}
LAMBDA_ONSET = {2: 392, 10: 392, 13: 391}

wb = openpyxl.load_workbook(SRC + "%R.xlsx", data_only=True)

def load_run(n):
    ws = wb[f"TiO2-{n}"]
    data = np.array([[c.value for c in row] for row in ws.iter_rows(min_row=2)], dtype=float)
    nm, pctR, Fr, hv, tauc_indirect = data[:, 0], data[:, 1], data[:, 3], data[:, 4], data[:, 7]
    order = np.argsort(nm)
    return nm[order], pctR[order], Fr[order], hv[order], tauc_indirect[order]

fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 3.1))

for n, meta in REPS.items():
    nm, pctR, Fr, hv, tauc = load_run(n)
    mask_vis = (nm >= 300) & (nm <= 800)
    axA.plot(nm[mask_vis], pctR[mask_vis], color=meta["color"], lw=1.3, label=meta["label"])

axA.set_xlabel("Wavelength (nm)")
axA.set_ylabel("Diffuse reflectance, %R")
axA.set_xlim(300, 800)
clean_axes(axA, grid_axis="y")
axA.legend(frameon=False, loc="lower right", fontsize=7.2)
axA.text(0.5, 1.04, "(a)", transform=axA.transAxes, fontsize=11, fontweight="bold", ha="center")

eg_text_lines = []
for n, meta in REPS.items():
    nm, pctR, Fr, hv, tauc = load_run(n)
    mask = (hv >= 2.6) & (hv <= 3.8)
    axB.plot(hv[mask], tauc[mask], color=meta["color"], lw=1.3, label=meta["label"])
    # local linear fit on the steepest-rising window just below the reported Eg
    eg = EG_INDIRECT[n]
    fit_mask = (hv >= eg - 0.12) & (hv <= eg + 0.03)
    if fit_mask.sum() >= 5:
        p = np.polyfit(hv[fit_mask], tauc[fit_mask], 1)
        x_ext = np.array([eg - 0.35, eg + 0.02])
        axB.plot(x_ext, np.polyval(p, x_ext), color=meta["color"], lw=0.9, ls="--", alpha=0.7)
    eg_text_lines.append((meta["label"].split(" (")[0], eg, meta["color"]))

axB.axhline(0, color=MUTED, lw=0.7)
axB.set_xlabel("Photon energy, h$\\nu$ (eV)")
axB.set_ylabel("[F(R)h$\\nu$]$^{1/2}$ (a.u.)")
axB.set_xlim(2.6, 3.8)
axB.set_ylim(bottom=-0.08)
clean_axes(axB, grid_axis="y")
axB.text(0.5, 1.04, "(b)", transform=axB.transAxes, fontsize=11, fontweight="bold", ha="center")
# compact Eg legend box (values are too close together to label on-axis)
box_y0 = 3.55
axB.text(2.68, box_y0 + 0.30, "$E_g$ (indirect):", fontsize=7.3, color=INK, ha="left")
for i, (name, eg, col) in enumerate(eg_text_lines):
    axB.text(2.68, box_y0 - i * 0.30, f"{name}: {eg:.2f} eV", fontsize=7.3, color=col, ha="left")

save(fig, "fig02_drs_tauc")

print("=" * 78)
print("DONE")
print("=" * 78)
