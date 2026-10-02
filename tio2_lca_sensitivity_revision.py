"""One-at-a-time sensitivity analysis of the single-process LCA (revision 2026-09-25).
The study uses laboratory grid electricity only (no photovoltaic supply), so the electricity source is not varied.

Baseline = tio2_lca_single_process.py (inventory per the Methods; FU = 1 g catalyst). Parameters varied one at a time:
  1. Drying-oven electricity  ±20 %  (9.600 kWh -> 7.680 / 11.520 kWh)
  2. Wastewater              total wastewater ±20 % (2.000 L -> 1.600 / 2.400 L); DI-water input (2.000 L) unchanged
  3. Catalyst yield           95 % instead of 100 % (every per-batch flow divided by 0.95); 95 % is the yield used in the
                              laboratory's earlier LCA workbook (TiO2-SolarModified_LCA_v4.xlsx, sheet Analysis, cell C1)
Writes ../data/lca_sensitivity_revision_results.json
"""
import copy, importlib.util, json, os

here = os.path.dirname(os.path.abspath(__file__))
os.chdir(here)
spec = importlib.util.spec_from_file_location('lca', os.path.join(here, 'tio2_lca_single_process.py'))
lca = importlib.util.module_from_spec(spec); spec.loader.exec_module(lca)      # re-writes the (identical) baseline json

CODES = lca.MIDPOINT_CODES
EXTRA = ['DALY', 'species_yr', 'USD2013', 'mPt_HH', 'mPt_ES', 'mPt_RS']

def evaluate(elec_kwh=None, wastewater_m3=None, yield_frac=1.0):
    fg = copy.deepcopy(lca.FOREGROUND)
    if elec_kwh is not None:
        fg['Electricity_grid_TH']['qty'] = elec_kwh
    if wastewater_m3 is not None:
        fg['Wastewater']['qty'] = wastewater_m3
    cf = copy.deepcopy(lca.CF)
    tot = {k: 0.0 for k in CODES + EXTRA}
    for flow in fg:
        amt = lca.FOREGROUND_KG.get(flow, fg[flow]['qty'])
        for k in CODES + EXTRA:
            tot[k] += amt * cf[flow][k]
    tot['GWP'] += lca.transport_packaging_gwp
    tot = {k: v / yield_frac for k, v in tot.items()}
    tot['single_score_mPt'] = tot['mPt_HH'] + tot['mPt_ES'] + tot['mPt_RS']
    tot['electricity_kWh'] = fg['Electricity_grid_TH']['qty'] / yield_frac
    return tot


E0 = lca.ELECTRICITY_TOTAL_KWH
DRY = next(v for k, v in lca.ENERGY_BREAKDOWN_KWH.items() if k.startswith('Drying'))
cases = [
    ('Baseline', 'Grid electricity; drying 9.600 kWh; wastewater 2.000 L; yield 100%', {}),
    ('Drying-oven electricity −20%', 'drying 7.680 kWh', dict(elec_kwh=E0 - 0.2 * DRY)),
    ('Drying-oven electricity +20%', 'drying 11.520 kWh', dict(elec_kwh=E0 + 0.2 * DRY)),
    ('Wastewater −20%', 'total wastewater 1.600 L', dict(wastewater_m3=0.001600)),
    ('Wastewater +20%', 'total wastewater 2.400 L', dict(wastewater_m3=0.002400)),
    ('Catalyst yield 95%', 'all flows ÷ 0.95', dict(yield_frac=0.95)),
]
res = []
base = None
for name, desc, kw in cases:
    t = evaluate(**kw)
    if base is None:
        base = t
    pct = {k: 100 * (t[k] / base[k] - 1) if base[k] else 0.0 for k in t}
    res.append({'case': name, 'setting': desc, 'totals': t, 'pct_change_vs_baseline': pct})
out = {'baseline_check_GWP': base['GWP'], 'cases': res}
json.dump(out, open('../data/lca_sensitivity_revision_results.json', 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
for r in res:
    t, p = r['totals'], r['pct_change_vs_baseline']
    print(f"{r['case']:34s} elec {t['electricity_kWh']:.3f} kWh | GWP {t['GWP']:.3f} ({p['GWP']:+.3f}%) | WC {t['WC']:.5f} ({p['WC']:+.3f}%) | "
          f"single score {t['single_score_mPt']:.3f} ({p['single_score_mPt']:+.3f}%) | max |Δ| over 18 midpoints {max(abs(p[c]) for c in CODES):.3f}%")
