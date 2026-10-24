"""Fig. 6 (formerly Fig. 4): ML benchmark, redrawn from manuscript_data.json with the styling of s05_fig03_05_rsm_surfaces.py
(FIGURE 2 block). Only change (author request 2026-09-24): panel (b) legend moved to the lower right, where it no
longer covers the negative RSM LOOCV bars."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

M = json.load(open('../data/manuscript_data.json', encoding='utf-8'))
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID, BASE = "#0b0b0b", "#52514e", "#e1e0d9", "#c3c2b7"
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "axes.edgecolor": BASE, "axes.linewidth": 0.8, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.labelcolor": INK, "figure.facecolor": "white",
    "savefig.facecolor": "white", "axes.axisbelow": True,
})


def clean_axes(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="y", color=GRID, linewidth=0.7, zorder=0)


ML = {r['Response']: r for r in M['ml_summary']}
RSM = {r['Response']: r for r in M['rsm_summary']}
model_names = ["GPR", "RF", "SVR", "KNR", "GBR", "MLP", "XGBR", "LGBMR", "DTR", "RSM"]
row = ML['Photodeg_pct']
insample = np.array([row[f"{n}_R2_insample"] for n in model_names])
loocv = np.array([row[f"{n}_R2_LOOCV"] for n in model_names])

fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.6, 3.1), gridspec_kw={"width_ratios": [1.55, 1.25]})
x = np.arange(len(model_names)); w = 0.38
axA.bar(x - w / 2, insample, w, color=BLUE, label="In-sample $R^2$", zorder=3)
axA.bar(x + w / 2, loocv, w, color=ORANGE, label="Leave-one-out CV $R^2$", zorder=3)
axA.axhline(0, color=INK2, linewidth=0.8, zorder=2)
axA.set_xticks(x); axA.set_xticklabels(model_names, rotation=35, ha="right")
for lbl in axA.get_xticklabels():
    if lbl.get_text() == "SVR":
        lbl.set_fontweight("bold")
axA.set_ylabel("$R^2$ (Photodegradation efficiency)")
axA.set_ylim(min(-1.6, float(loocv.min()) - 0.15), 1.08)
clean_axes(axA); axA.legend(frameon=False, loc="lower left", fontsize=7.3)
axA.set_title("(a)", loc="center", fontsize=9, fontweight="bold")

resps = ["Photodeg_pct", "k_app_min-1", "AQE_pct", "EEO_kWh_m3_order", "SEC_total_Wh_mg"]
resp_labels = ["Photodeg\n(%)", "$k_{app}$\n(min$^{-1}$)", "AQE\n(%)", "EEO\n(kWh m$^{-3}$\nord$^{-1}$)", "SEC$_{total}$\n(Wh mg$^{-1}$)"]
rsm_in = [RSM[r]["R2"] for r in resps]
rsm_loo = [ML[r]["RSM_R2_LOOCV"] for r in resps]
svr_loo = [ML[r]["SVR_R2_LOOCV"] for r in resps]
xr = np.arange(len(resps)); w2 = 0.27
axB.bar(xr - w2, rsm_in, w2, color=BLUE, label="RSM in-sample", zorder=3)
axB.bar(xr, rsm_loo, w2, color=ORANGE, label="RSM LOOCV", zorder=3)
axB.bar(xr + w2, svr_loo, w2, color=AQUA, label="SVR LOOCV (best)", zorder=3)
axB.axhline(0, color=INK2, linewidth=0.8, zorder=2)
axB.set_xticks(xr); axB.set_xticklabels(resp_labels, fontsize=6.6)
axB.set_ylabel("$R^2$")
axB.set_ylim(-1.3, 1.08)
clean_axes(axB); axB.legend(frameon=False, fontsize=6.4, loc="lower right")
axB.set_title("(b)", loc="center", fontsize=9, fontweight="bold")
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"../figures/fig06_ml_benchmark.{ext}", dpi=600, bbox_inches="tight")
print("saved fig06_ml_benchmark")
