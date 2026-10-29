"""Lack-of-fit test for the photodegradation RSM model and design-support check (review revision 2026-10-24).

1. Refits the quadratic model for photodegradation (coded units, ordinary least squares) to the 15 design runs and tests
   lack of fit against pure error from the three design-center replicates (Section 3.2).
2. Builds the convex hull of the 13 distinct tested design points. Run 4 was performed at the face center (+1, 0, 0)
   instead of the Box-Behnken edge point (+1, +1, 0), so the hull is the full-BBD hull cut by x1 + x2 <= 1
   (Sections 2.4, 2.11). Reports x1 + x2 and hull membership for the validated optimum and V4-V7 (Sections 3.4, 3.5).

Reads ../data/manuscript_data.json (written by s02); writes ../data/design_support_results.json.
"""
import json

import numpy as np
from scipy import stats
from scipy.spatial import ConvexHull, Delaunay

M = json.load(open('../data/manuscript_data.json', encoding='utf-8'))
bbd = M['bbd']


def coded(g, c, t):
    return np.array([(g - 1.25) / 0.75, (c - 0.05) / 0.03, (t - 150) / 30])


X = np.array([coded(b['Glucose_g'], b['CuSO4_g'], b['Time_min']) for b in bbd])
y = np.array([b['Photodeg_pct'] for b in bbd])
x1, x2, x3 = X.T
F = np.column_stack([np.ones(len(y)), x1, x2, x3, x1**2, x2**2, x3**2, x1*x2, x1*x3, x2*x3])
beta, *_ = np.linalg.lstsq(F, y, rcond=None)
sse = float(((y - F @ beta) ** 2).sum())
sst = float(((y - y.mean()) ** 2).sum())

center = np.all(np.isclose(X, 0), axis=1)
pe = float(((y[center] - y[center].mean()) ** 2).sum())
df_pe = int(center.sum()) - 1
df_lof = len(y) - F.shape[1] - df_pe
F_lof = ((sse - pe) / df_lof) / (pe / df_pe)
p_lof = float(stats.f.sf(F_lof, df_lof, df_pe))
print(f'R2 = {1 - sse / sst:.3f}; lack of fit F = {F_lof:.1f} (df {df_lof}, {df_pe}), p = {p_lof:.4f}; '
      f'pure-error SD = {np.sqrt(pe / df_pe):.3f} percentage points')

pts = np.unique(X.round(9), axis=0)
hull = ConvexHull(pts)
tri = Delaunay(pts)
facets = np.unique(hull.equations.round(6), axis=0)
points = [('Validated optimum', M['new_optimum_condition'])] + \
         [(v['Set'], v) for v in M['val_summary'] if v['Set'] in ('V4', 'V5', 'V6', 'V7')]
support = []
for name, c in points:
    x = coded(c['Glucose_g'], c['CuSO4_g'], c['Time_min'])
    row = dict(point=name, Glucose_g=c['Glucose_g'], CuSO4_g=c['CuSO4_g'], Time_min=c['Time_min'],
               coded=[float(v) for v in x], L1=float(np.abs(x).sum()), x1_plus_x2=float(x[0] + x[1]),
               inside_tested_hull=bool(tri.find_simplex(x) >= 0))
    support.append(row)
    print(f"{name:18s} x1+x2 = {row['x1_plus_x2']:.3f}  L1 = {row['L1']:.3f}  inside tested hull: {row['inside_tested_hull']}")
# CuSO4 at which the tested hull ends for the validated glucose dose and time (x1 + x2 = 1)
x_opt = coded(1.57, 0.072, 140)
cu_edge = 0.05 + 0.03 * (1 - x_opt[0])
print(f'tested hull boundary at 1.570 g glucose, 140 min: CuSO4 = {cu_edge:.4f} g')

out = {'lack_of_fit': {'R2': 1 - sse / sst, 'SSE': sse, 'pure_error_SS': pe, 'df_lack_of_fit': df_lof, 'df_pure_error': df_pe,
                       'F': F_lof, 'p': p_lof, 'pure_error_sd_pct_points': float(np.sqrt(pe / df_pe)),
                       'coefficients': [float(b) for b in beta]},
       'tested_hull': {'n_distinct_points': int(len(pts)), 'volume': float(hull.volume),
                       'facets_n_dot_x_plus_b_le_0': facets.tolist(),
                       'extra_constraint_vs_full_BBD': 'x1 + x2 <= 1',
                       'CuSO4_boundary_at_validated_glucose_time_g': cu_edge},
       'points': support}
json.dump(out, open('../data/design_support_results.json', 'w', encoding='utf-8'), indent=1)
print('wrote ../data/design_support_results.json')
