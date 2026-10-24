"""Fig. 7 (formerly Fig. 5): trade-off between photodegradation and EEO / SEC_total, redrawn (revision 2026-09-24)
with the NSGA-II front constrained to the convex hull of the Box–Behnken design points (s07_optimization_tested_hull.py).
The desirability, weighted-sum and Tchebycheff scalarizations saturate (D = 1) over a region, so that region is
drawn (predicted responses of points sampled inside it) instead of three single markers. Styling as s05_fig03_05_rsm_surfaces.py.
Review revision 2026-10-24: the region and front are restricted to the convex hull of the tested points (BBD hull and
x1 + x2 <= 1, since Run 4 was performed at (+1, 0, 0))."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

M = json.load(open('../data/manuscript_data.json', encoding='utf-8'))
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
INK, INK2, MUTED, GRID, BASE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "axes.edgecolor": BASE, "axes.linewidth": 0.8, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.labelcolor": INK, "figure.facecolor": "white",
    "savefig.facecolor": "white", "axes.axisbelow": True,
})
TERMS = ['Intercept', 'X1', 'X2', 'X3', 'I(X1 ** 2)', 'I(X2 ** 2)', 'I(X3 ** 2)', 'X1:X2', 'X1:X3', 'X2:X3']
RESP = ['Photodeg_pct', 'EEO_kWh_m3_order', 'SEC_total_Wh_mg']
COEF = {r: np.array([M['rsm_coefs'][r][t] for t in TERMS]) for r in RESP}


def predict(r, x1, x2, x3):
    return np.column_stack([np.ones_like(x1), x1, x2, x3, x1**2, x2**2, x3**2, x1*x2, x1*x3, x2*x3]) @ COEF[r]


bbd = M['bbd']
best = {'Photodeg_pct': max(b['Photodeg_pct'] for b in bbd), 'EEO_kWh_m3_order': min(b['EEO_kWh_m3_order'] for b in bbd),
        'SEC_total_Wh_mg': min(b['SEC_total_Wh_mg'] for b in bbd)}
rng = np.random.default_rng(0)
X = rng.uniform(-1, 1, (400000, 3)); X = X[(np.abs(X).sum(1) <= 2) & (X[:, 0] + X[:, 1] <= 1)]
pr = {r: predict(r, *X.T) for r in RESP}
plat = (pr['Photodeg_pct'] >= best['Photodeg_pct']) & (pr['EEO_kWh_m3_order'] <= best['EEO_kWh_m3_order']) & (pr['SEC_total_Wh_mg'] <= best['SEC_total_Wh_mg'])
idx = np.flatnonzero(plat); idx = idx[rng.choice(len(idx), min(1500, len(idx)), replace=False)]
center = [b for b in bbd if b['Run'] in (13, 14, 15)]
pareto = M['pareto']; nop = M['new_optimum_model_prediction']

fig, (axL, axR) = plt.subplots(1, 2, figsize=(9.0, 3.9))
for ax, col, lab, ylim in ((axL, 'EEO_kWh_m3_order', "EEO (kWh m$^{-3}$ order$^{-1}$)", (-3000, 62000)),
                           (axR, 'SEC_total_Wh_mg', "SEC$_{total}$ (Wh mg$^{-1}$)", (-1000, 38000))):
    ax.scatter([b['Photodeg_pct'] for b in bbd], [b[col] for b in bbd], s=30, color=MUTED, edgecolor=INK2, linewidth=0.5,
               label=f"Design points (n = {len(bbd)})", zorder=2)
    ax.scatter(pr['Photodeg_pct'][idx], pr[col][idx], s=4, color=ORANGE, alpha=0.25, linewidth=0, zorder=3,
               label="Scalarization plateau (D = 1)")
    ax.scatter([p['Photodeg_pct'] for p in pareto], [p[col] for p in pareto], s=12, color=BLUE, alpha=0.7,
               label="NSGA-II Pareto front", zorder=4)
    ax.scatter(M['center_photodeg_mean'], np.mean([b[col] for b in center]), marker="*", s=200, color=AQUA, edgecolor=INK,
               linewidth=0.9, zorder=6, label="Design center\n(Runs 13–15, n = 3)")
    ax.scatter(nop['Photodeg_pct'], nop[col], marker="P", s=160, color=VIOLET, edgecolor=INK, linewidth=0.9, zorder=6,
               label="New optimum\n(experimentally validated, n = 9)")
    ax.set_xlabel("Photodegradation efficiency (%)"); ax.set_ylabel(lab); ax.set_ylim(*ylim)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="y", color=GRID, linewidth=0.7, zorder=0)
axL.set_title("(a)", loc="center", fontsize=9, fontweight="bold")
axR.set_title("(b)", loc="center", fontsize=9, fontweight="bold")
leg = axR.legend(frameon=True, facecolor="white", edgecolor=BASE, framealpha=0.95, fontsize=6.4, loc="upper right", markerscale=1.0)
for h in leg.legend_handles:
    try:
        h.set_alpha(1)
    except Exception:
        pass
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"../figures/fig07_pareto_optimization.{ext}", dpi=600, bbox_inches="tight")
print("saved fig07_pareto_optimization; plateau points drawn:", len(idx))
