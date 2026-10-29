"""Multi-objective optimization re-run under the Box–Behnken convex-hull constraint (revision 2026-09-24).

Reproduces Section 5 of s02_rsm_ml_validation_kinetics.py exactly (same fitted RSM coefficients, desirability scaling, differential
evolution and NSGA-II settings, seed 42), changing only the design-support constraint:

  before:  ||x||_2 <= sqrt(2)                      (sphere through the design points)
  now:     |x1| + |x2| + |x3| <= 2,  |xi| <= 1     (convex hull of the 3-factor BBD = cuboctahedron)

The sphere contains regions outside the convex hull of the tested points (e.g. the old NSGA-II best-EEO solution,
L1 = 2.124), so it does not guarantee interpolation. Run with  --check  to reproduce the old spherical results first.
Reads ../data/manuscript_data.json; writes ../data/optimization_hull_results.json and updates opt_summary / pareto in
../data/manuscript_data.json (the spherical results are kept under opt_summary_sphere / pareto_sphere).

Review revision 2026-10-24: Run 4 was performed at the face centre (+1, 0, 0) instead of the BBD edge point (+1, +1, 0),
so the convex hull of the points actually tested is the BBD hull cut by  x1 + x2 <= 1.  That is now the default
constraint (MODE 'tested'); run with  --bbd  to reproduce the full-BBD hull results of revision 2026-09-24.
Writes ../data/optimization_tested_hull_results.json; the full-BBD results are kept in manuscript_data.json under
opt_summary_bbd_hull / pareto_bbd_hull.
"""
import json, sys
import numpy as np
from scipy.optimize import differential_evolution

DATA = '../data/manuscript_data.json'
M = json.load(open(DATA, encoding='utf-8'))
RNG = 42
BOUNDS_REAL = [(0.5, 2.0), (0.02, 0.08), (120, 180)]
RESP = ['Photodeg_pct', 'EEO_kWh_m3_order', 'SEC_total_Wh_mg']
TERMS = ['Intercept', 'X1', 'X2', 'X3', 'I(X1 ** 2)', 'I(X2 ** 2)', 'I(X3 ** 2)', 'X1:X2', 'X1:X3', 'X2:X3']
COEF = {r: np.array([M['rsm_coefs'][r][t] for t in TERMS]) for r in RESP}
bbd = M['bbd']
lo = {r: min(b[r] for b in bbd) for r in RESP}
hi = {r: max(b[r] for b in bbd) for r in RESP}


def code(x, a, b):
    return (np.asarray(x, dtype=float) - (a + b) / 2) / ((b - a) / 2)


def coded(g, c, t):
    return code(g, 0.5, 2.0), code(c, 0.02, 0.08), code(t, 120, 180)


def predict(resp, g, c, t):
    x1, x2, x3 = (np.atleast_1d(v) for v in coded(g, c, t))
    F = np.column_stack([np.ones_like(x1), x1, x2, x3, x1**2, x2**2, x3**2, x1*x2, x1*x3, x2*x3])
    p = F @ COEF[resp]
    return float(p[0]) if np.ndim(g) == 0 else p


MODE = 'sphere' if '--check' in sys.argv else 'hull' if '--bbd' in sys.argv else 'tested'


def excess(g, c, t):
    x = np.array(coded(g, c, t), dtype=float)
    if MODE == 'sphere':
        return np.sqrt((x**2).sum(axis=0)) - np.sqrt(2.0)
    bbd_hull = np.abs(x).sum(axis=0) - 2.0      # |xi| <= 1 is enforced by the variable bounds
    if MODE == 'hull':
        return bbd_hull
    return np.maximum(bbd_hull, x[0] + x[1] - 1.0)   # Run 4 at (+1, 0, 0): the (+1, +1, 0) edge was not tested


def penalty(g, c, t, weight=50.0):
    return weight * max(0.0, float(excess(g, c, t)))**2


def desir(g, c, t):
    d = predict('Photodeg_pct', g, c, t); e = predict('EEO_kWh_m3_order', g, c, t); s = predict('SEC_total_Wh_mg', g, c, t)
    return (np.clip((d - lo['Photodeg_pct']) / (hi['Photodeg_pct'] - lo['Photodeg_pct']), 1e-6, 1),
            np.clip((hi['EEO_kWh_m3_order'] - e) / (hi['EEO_kWh_m3_order'] - lo['EEO_kWh_m3_order']), 1e-6, 1),
            np.clip((hi['SEC_total_Wh_mg'] - s) / (hi['SEC_total_Wh_mg'] - lo['SEC_total_Wh_mg']), 1e-6, 1))


def f_des(p):
    a, b, c = desir(*p); return -(a * b * c) ** (1 / 3) + penalty(*p)


def f_ws(p):
    a, b, c = desir(*p); return -(a + b + c) / 3 + penalty(*p)


def f_tch(p, ideal):
    a, b, c = desir(*p); return max(abs(ideal[0] - a), abs(ideal[1] - b), abs(ideal[2] - c)) + penalty(*p)


kw = dict(seed=RNG, popsize=15, maxiter=300, tol=1e-10)
sol = {'Desirability (DE)': differential_evolution(f_des, BOUNDS_REAL, **kw).x,
       'Weighted sum (equal)': differential_evolution(f_ws, BOUNDS_REAL, **kw).x,
       'Tchebycheff': differential_evolution(f_tch, BOUNDS_REAL, args=((1, 1, 1),), **kw).x}

from pymoo.core.problem import Problem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize as pymoo_minimize


class P3(Problem):
    def __init__(self):
        super().__init__(n_var=3, n_obj=3, n_ieq_constr=1, xl=np.array([0.5, 0.02, 120]), xu=np.array([2.0, 0.08, 180]))

    def _evaluate(self, X, out, *a, **k):
        g, c, t = X[:, 0], X[:, 1], X[:, 2]
        out['F'] = np.column_stack([-predict('Photodeg_pct', g, c, t), predict('EEO_kWh_m3_order', g, c, t), predict('SEC_total_Wh_mg', g, c, t)])
        out['G'] = np.column_stack([excess(g, c, t)])


res = pymoo_minimize(P3(), NSGA2(pop_size=100), ('n_gen', 150), seed=RNG, verbose=False)
pareto = [dict(Glucose_g=float(x[0]), CuSO4_g=float(x[1]), Time_min=float(x[2]),
               Photodeg_pct=float(-f[0]), EEO_kWh_m3_order=float(f[1]), SEC_total_Wh_mg=float(f[2])) for x, f in zip(res.X, res.F)]
best_deg = max(pareto, key=lambda r: r['Photodeg_pct'])
best_eeo = min(pareto, key=lambda r: r['EEO_kWh_m3_order'])


def row(name, g, c, t):
    return dict(Method=name, Glucose_g=float(g), CuSO4_g=float(c), Time_min=float(t),
                Photodeg_pct=predict('Photodeg_pct', g, c, t), EEO_kWh_m3_order=predict('EEO_kWh_m3_order', g, c, t),
                SEC_total_Wh_mg=predict('SEC_total_Wh_mg', g, c, t))


opt = [row(k, *v) for k, v in sol.items()] + [dict(Method='NSGA-II (best degradation)', **best_deg), dict(Method='NSGA-II (best EEO)', **best_eeo)]
L1 = lambda r: float(np.abs(np.array(coded(r['Glucose_g'], r['CuSO4_g'], r['Time_min']))).sum())
L2 = lambda r: float(np.sqrt((np.array(coded(r['Glucose_g'], r['CuSO4_g'], r['Time_min']))**2).sum()))
print(f'mode = {MODE}; NSGA-II front: {len(pareto)} solutions; max L1 on front = {max(L1(p) for p in pareto):.4f}; max L2 = {max(L2(p) for p in pareto):.4f}')
for r in opt:
    print(f"  {r['Method']:28s} {r['Glucose_g']:.4f} {r['CuSO4_g']:.5f} {r['Time_min']:.3f} | deg {r['Photodeg_pct']:.3f} EEO {r['EEO_kWh_m3_order']:.1f} SEC {r['SEC_total_Wh_mg']:.1f} | L1 {L1(r):.3f} L2 {L2(r):.3f}")

if MODE == 'sphere':
    old = {r['Method']: r for r in M['opt_summary']}
    for r in opt:
        o = old[r['Method']]
        print('  reproduces stored', r['Method'], all(abs(r[k] - o[k]) < 1e-3 * max(1, abs(o[k])) for k in ('Glucose_g', 'CuSO4_g', 'Time_min', 'Photodeg_pct')))
    sys.exit(0)

CONSTRAINT = {'hull': '|x1|+|x2|+|x3| <= 2 and |xi| <= 1 (convex hull of the 3-factor Box-Behnken design)',
              'tested': '|x1|+|x2|+|x3| <= 2, x1 + x2 <= 1 and |xi| <= 1 (convex hull of the 13 distinct tested design points)'}[MODE]
out = {'constraint': CONSTRAINT,
       'opt_summary': opt, 'pareto': pareto,
       'validated_point_L1': L1(dict(Glucose_g=1.57, CuSO4_g=0.072, Time_min=140)),
       'validated_point_x1_plus_x2': float(sum(coded(1.57, 0.072, 140)[:2]))}


# ---- D = 1 region (Table 5 ranges): points of the hull where every desirability saturates, i.e. predicted
# photodegradation >= best observed run and predicted EEO, SEC_total <= best observed runs.
def desirability_plateau():
    from scipy.optimize import minimize
    n = 161                                                     # grid screen of the factor box
    g, c, t = (a.ravel() for a in np.meshgrid(np.linspace(0.5, 2.0, n), np.linspace(0.02, 0.08, n), np.linspace(120, 180, n), indexing='ij'))
    xg = np.array(coded(g, c, t))
    hull = np.abs(xg).sum(0) <= 2 + 1e-12
    if MODE == 'tested':
        hull &= xg[0] + xg[1] <= 1 + 1e-12
    d, e, s = predict('Photodeg_pct', g, c, t), predict('EEO_kWh_m3_order', g, c, t), predict('SEC_total_Wh_mg', g, c, t)
    plat = hull & (d >= hi['Photodeg_pct']) & (e <= lo['EEO_kWh_m3_order']) & (s <= lo['SEC_total_Wh_mg'])
    mm = lambda v: [float(v[plat].min()), float(v[plat].max())]
    res = {'grid_points_per_axis': n, 'hull_fraction_of_box': float(hull.mean()), 'plateau_fraction_of_hull': float(plat[hull].mean()),
           'glucose_g': mm(g), 'CuSO4_g': mm(c), 'time_min': mm(t), 'pred_photodeg_pct': mm(d), 'pred_EEO': mm(e), 'pred_SEC': mm(s),
           'observed_best': {'photodeg': hi['Photodeg_pct'], 'EEO': lo['EEO_kWh_m3_order'], 'SEC': lo['SEC_total_Wh_mg']},
           'min_pred_in_hull': {'EEO': float(e[hull].min()), 'SEC': float(s[hull].min())},
           'min_pred_in_box': {'EEO': float(e.min()), 'SEC': float(s.min())},
           'negative_outside_hull_points': int(((e < 0) | (s < 0))[~hull].sum())}
    # exact extent (reported in Table 5): SLSQP from 40 feasible grid starts, coded space
    dec = lambda x: (1.25 + 0.75 * x[0], 0.05 + 0.03 * x[1], 150 + 30 * x[2])
    F = [lambda x: predict('Photodeg_pct', *dec(x)) - hi['Photodeg_pct'],
         lambda x: lo['EEO_kWh_m3_order'] - predict('EEO_kWh_m3_order', *dec(x)),
         lambda x: lo['SEC_total_Wh_mg'] - predict('SEC_total_Wh_mg', *dec(x))]
    for k in range(8):                                          # |x1|+|x2|+|x3| <= 2 as eight linear inequalities
        sg = np.array([1 if k & 1 else -1, 1 if k & 2 else -1, 1 if k & 4 else -1])
        F.append(lambda x, sg=sg: 2 - sg @ x)
    if MODE == 'tested':
        F.append(lambda x: 1 - x[0] - x[1])
    cons = [{'type': 'ineq', 'fun': f} for f in F]
    feas = lambda x: all(float(np.squeeze(f(x))) >= -1e-7 for f in F) and np.all(np.abs(x) <= 1 + 1e-9)
    G = np.array(np.meshgrid(*[np.linspace(-1, 1, 41)] * 3, indexing='ij')).reshape(3, -1).T
    starts = [x for x in G if feas(x)]
    starts = [starts[i] for i in np.random.default_rng(0).choice(len(starts), min(40, len(starts)), replace=False)]

    def extreme(obj):
        out_ = []
        for sign in (1, -1):
            vals = []
            for x0 in starts:
                r = minimize(lambda x: sign * float(np.squeeze(obj(x))), x0, constraints=cons, bounds=[(-1, 1)] * 3,
                             method='SLSQP', options={'ftol': 1e-12, 'maxiter': 1000})
                if feas(r.x):
                    vals.append(float(np.squeeze(obj(r.x))))
            out_.append(min(vals) if sign == 1 else max(vals))
        return out_
    ex = {'glucose_g': extreme(lambda x: dec(x)[0]), 'CuSO4_g': extreme(lambda x: dec(x)[1]), 'time_min': extreme(lambda x: dec(x)[2])}
    for resp in RESP:
        ex['pred_' + resp] = extreme(lambda x, r=resp: predict(r, *dec(x)))
    res['exact_bounds_SLSQP'] = ex
    return res


out['desirability_plateau'] = desirability_plateau()
pb = out['desirability_plateau']['exact_bounds_SLSQP']
print('D = 1 region:', {k: [round(v, 3) for v in pb[k]] for k in pb})
out_file = {'hull': 'optimization_hull_results.json', 'tested': 'optimization_tested_hull_results.json'}[MODE]
json.dump(out, open('../data/' + out_file, 'w', encoding='utf-8'), indent=1)
if 'opt_summary_sphere' not in M:
    M['opt_summary_sphere'], M['pareto_sphere'] = M['opt_summary'], M['pareto']
if MODE == 'tested':                            # full-BBD results from the preceding  --bbd  run (run_all.py step 7a)
    B = json.load(open('../data/optimization_hull_results.json', encoding='utf-8'))
    M['opt_summary_bbd_hull'], M['pareto_bbd_hull'] = B['opt_summary'], B['pareto']
M['opt_summary'], M['pareto'] = opt, pareto
M['optimization_constraint'] = out['constraint']
json.dump(M, open(DATA, 'w', encoding='utf-8'), indent=1)
print(f'wrote ../data/{out_file} and updated ../data/manuscript_data.json')
