"""Uncertainty ranges and process-improvement scenarios for the single-process LCA (revision 2026-10-25, ET&I).

Functional unit: one laboratory batch initiated with 1.000 g of commercial anatase TiO2 (s11_lca_baseline.py). The
baseline electricity inventory is a rated-power upper bound (equipment-rated power x operating time).

Modelled scenarios, not new experiments: each one changes only the electricity inventory and keeps every other flow,
the system boundary and the 100 % recovery assumption unchanged.
  S1  Idealized full-load drying allocation: one 12 h, 800 W drying cycle serves N equivalent one-gram batches
      (N = 5, 10, 20), so the drying electricity per batch is 9.600 / N kWh. Upper-bound reduction: it assumes that oven
      power and drying time do not change with load.
  S2  Drying-cycle redesign: drying electricity -25 % and -50 % (shorter cycle or lower set-point), assuming the
      solid still reaches constant mass; not tested experimentally.
  S3  Electricity-supply substitution: 50 % and 100 % of the process electricity from rooftop photovoltaics.
      Characterization factors for 'Electricity, low voltage {TH} | electricity production, photovoltaic, 3 kWp' are
      taken from the laboratory LCA workbook (TiO2-SolarModified_LCA_v4.xlsx, sheet CF, row 7), rounded to three
      significant figures like the grid factors in s11_lca_baseline.py. No storage is included; temporal matching or
      grid-mediated accounting is assumed.
  S4  Combined: S1 (N = 10) + S2 (-50 %) + S3 (100 % PV).
Uncertainty ranges (one parameter at a time, then two bounding combinations):
  drying duty factor 0.25-1.00 (E_dry = P_rated x t x f_duty; 1.00 = rated-power upper bound), recovery yield 80-100 %
  (every flow divided by the yield), drying allocation 1-10 batches, electricity characterization factors +/-20 %
  (all categories), electricity of the non-drying equipment +/-20 %.
Also writes the illustrative production GWP per mg of methylene blue photodegraded in one test (50 mg catalyst,
0.100 L of 2.500 mg/L MB, 22.304 % photodegradation) for 1, 5 and 10 hypothetical uses, and the number of uses at which
the production GWP allocated to one test equals the GWP of the test's own LED electricity (100 W, 60 min, grid).
Writes ../data/lca_scenario_results.json
"""
import importlib.util, json, os

here = os.path.dirname(os.path.abspath(__file__))
os.chdir(here)
spec = importlib.util.spec_from_file_location('lca', os.path.join(here, 's11_lca_baseline.py'))
lca = importlib.util.module_from_spec(spec); spec.loader.exec_module(lca)      # re-writes the (identical) baseline json

CODES = lca.MIDPOINT_CODES
EXTRA = ['DALY', 'species_yr', 'USD2013', 'mPt_HH', 'mPt_ES', 'mPt_RS']
CF_PV = {  # TiO2-SolarModified_LCA_v4.xlsx, sheet CF, row 7: low-voltage electricity, Thailand, rooftop PV 3 kWp
    "GWP": 0.0645, "SOD": 2.58e-08, "IR": 0.00473, "OF-HH": 0.000181, "FPMF": 0.000155, "OF-TE": 0.000191,
    "TA": 0.000347, "FE": 4.4e-05, "ME": 3.35e-06, "TET": 2.27, "FET": 0.0303, "MET": 0.0385, "HCT": 0.00922,
    "HNT": 0.272, "LU": 0.00195, "MRS": 0.00108, "FRS": 0.0163, "WC": 0.00226,
    "DALY": 2.54e-07, "species_yr": 4.03e-10, "USD2013": 0.0041, "mPt_HH": 4.24, "mPt_ES": 0.109, "mPt_RS": 0.0292,
}
DRY_KEY = next(k for k in lca.ENERGY_BREAKDOWN_KWH if k.startswith('Drying'))
DRY0 = lca.ENERGY_BREAKDOWN_KWH[DRY_KEY]


def evaluate(drying_kwh=DRY0, pv_share=0.0, other_scale=1.0, cf_scale=1.0, yield_frac=1.0):
    ops = {k: (drying_kwh if k == DRY_KEY else v * other_scale) for k, v in lca.ENERGY_BREAKDOWN_KWH.items()}
    elec = sum(ops.values())
    tot = {k: 0.0 for k in CODES + EXTRA}
    elec_part = {k: 0.0 for k in CODES + EXTRA}
    for flow in lca.FOREGROUND:
        if flow == 'Electricity_grid_TH':
            for k in CODES + EXTRA:
                v = elec * cf_scale * ((1 - pv_share) * lca.CF[flow][k] + pv_share * CF_PV[k])
                tot[k] += v
                elec_part[k] += v
            continue
        amt = lca.FOREGROUND_KG.get(flow, lca.FOREGROUND[flow]['qty'])
        for k in CODES + EXTRA:
            tot[k] += amt * lca.CF[flow][k]
    tot['GWP'] += lca.transport_packaging_gwp
    tot = {k: v / yield_frac for k, v in tot.items()}
    elec_part = {k: v / yield_frac for k, v in elec_part.items()}
    for d in (tot, elec_part):
        d['single_score_mPt'] = d['mPt_HH'] + d['mPt_ES'] + d['mPt_RS']
    return {'electricity_kWh': elec / yield_frac, 'drying_kWh': drying_kwh / yield_frac, 'pv_share': pv_share,
            'totals': tot, 'electricity_contribution': elec_part}


def with_change(r, base):
    r['pct_change_vs_baseline'] = {k: 100 * (r['totals'][k] / base['totals'][k] - 1) if base['totals'][k] else 0.0
                                   for k in r['totals']}
    r['electricity_share_of_GWP_pct'] = 100 * r['electricity_contribution']['GWP'] / r['totals']['GWP']
    return r


base = evaluate()
cases = [
    ('S0', 'Baseline (rated-power upper bound)', 'grid electricity; one batch per 12 h drying cycle', {}),
    ('S1a', 'Idealized full-load allocation, 5 batches', 'drying 9.600/5 kWh per batch', dict(drying_kwh=DRY0 / 5)),
    ('S1b', 'Idealized full-load allocation, 10 batches', 'drying 9.600/10 kWh per batch', dict(drying_kwh=DRY0 / 10)),
    ('S1c', 'Idealized full-load allocation, 20 batches', 'drying 9.600/20 kWh per batch', dict(drying_kwh=DRY0 / 20)),
    ('S2a', 'Drying-cycle redesign, −25%', 'drying 7.200 kWh per batch', dict(drying_kwh=DRY0 * 0.75)),
    ('S2b', 'Drying-cycle redesign, −50%', 'drying 4.800 kWh per batch', dict(drying_kwh=DRY0 * 0.50)),
    ('S3a', 'Electricity substitution, 50% rooftop PV', 'half of the process electricity from rooftop PV', dict(pv_share=0.5)),
    ('S3b', 'Electricity substitution, 100% rooftop PV', 'all process electricity from rooftop PV', dict(pv_share=1.0)),
    ('S4', 'Combined (S1b + S2b + S3b)', '10 batches per cycle, drying −50%, 100% rooftop PV',
     dict(drying_kwh=DRY0 * 0.50 / 10, pv_share=1.0)),
]
res = [dict(id=sid, scenario=name, setting=desc, **with_change(evaluate(**kw), base)) for sid, name, desc, kw in cases]

# uncertainty ranges: one parameter at a time around the rated-power upper bound, then two bounding combinations
ranges = []
for name, low, high in [
        ('Drying duty factor (0.25–1.00)', dict(drying_kwh=DRY0 * 0.25), dict(drying_kwh=DRY0 * 1.00)),
        ('Recovery yield (100–80%)', dict(yield_frac=1.0), dict(yield_frac=0.80)),
        ('Drying allocation (10–1 batches)', dict(drying_kwh=DRY0 / 10), dict(drying_kwh=DRY0)),
        ('Electricity characterization factors (−20% to +20%)', dict(cf_scale=0.8), dict(cf_scale=1.2)),
        ('Non-drying equipment electricity (−20% to +20%)', dict(other_scale=0.8), dict(other_scale=1.2)),
        ('Bounding combination', dict(drying_kwh=DRY0 * 0.25 / 10, yield_frac=1.0, cf_scale=0.8, other_scale=0.8),
         dict(drying_kwh=DRY0, yield_frac=0.80, cf_scale=1.2, other_scale=1.2))]:
    lo, hi = evaluate(**low), evaluate(**high)
    ranges.append({'parameter': name, 'low': {'electricity_kWh': lo['electricity_kWh'], 'GWP': lo['totals']['GWP'],
                                              'single_score_mPt': lo['totals']['single_score_mPt']},
                   'high': {'electricity_kWh': hi['electricity_kWh'], 'GWP': hi['totals']['GWP'],
                            'single_score_mPt': hi['totals']['single_score_mPt']}})
duty = [{'duty_factor': f, 'drying_kWh': DRY0 * f, **{k: evaluate(drying_kwh=DRY0 * f)['totals'][k] for k in ('GWP', 'single_score_mPt')},
         'electricity_kWh': evaluate(drying_kwh=DRY0 * f)['electricity_kWh']} for f in (0.25, 0.50, 0.75, 1.00)]

# functional linkage (production stage only; reuse hypothetical)
CAT_G, BATCH_G, MB_MG = 0.050, 1.000, 2.5 * 0.100 * 0.22304
TEST_KWH = 0.100 * 1.0                                        # 100 W LED, 60 min
test_gwp = TEST_KWH * lca.CF['Electricity_grid_TH']['GWP']
per_mg = {sid: {n: r['totals']['GWP'] * CAT_G / BATCH_G / n / MB_MG for n in (1, 5, 10)}
          for sid, r in ((x['id'], x) for x in res) if sid in ('S0', 'S4')}
break_even = base['totals']['GWP'] * CAT_G / BATCH_G / test_gwp
out = {'functional_unit': 'one laboratory batch initiated with 1.000 g of commercial anatase TiO2',
       'note': 'modelled scenarios; every flow except electricity is unchanged from the baseline',
       'pv_characterization_factors': CF_PV, 'baseline_check_GWP': base['totals']['GWP'],
       'scenarios': res, 'uncertainty_ranges': ranges, 'drying_duty_factor': duty,
       'illustrative_GWP_per_mg_MB_photodegraded': {'catalyst_g_per_test': CAT_G, 'MB_mg_photodegraded_per_test': MB_MG,
                                                     'kgCO2e_per_mg_by_scenario_and_uses': per_mg},
       'test_stage': {'LED_kWh_per_test': TEST_KWH, 'GWP_kgCO2e_per_test': test_gwp,
                      'production_GWP_per_test_single_use': base['totals']['GWP'] * CAT_G / BATCH_G,
                      'uses_for_production_GWP_to_equal_test_GWP': break_even}}
json.dump(out, open('../data/lca_scenario_results.json', 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
for r in res:
    t, p = r['totals'], r['pct_change_vs_baseline']
    print(f"{r['id']:4s} {r['scenario']:44s} elec {r['electricity_kWh']:6.3f} | GWP {t['GWP']:.3f} ({p['GWP']:+.3f}%) | "
          f"SS {t['single_score_mPt']:.3f} ({p['single_score_mPt']:+.3f}%) | IR {p['IR']:+.3f}% | TET {p['TET']:+.3f}%")
for r in ranges:
    print(f"{r['parameter']:52s} GWP {r['low']['GWP']:.3f}–{r['high']['GWP']:.3f} | SS {r['low']['single_score_mPt']:.3f}–{r['high']['single_score_mPt']:.3f}")
print('duty:', [(d['duty_factor'], round(d['electricity_kWh'], 3), round(d['GWP'], 3)) for d in duty])
print('per mg:', {k: {n: round(v, 3) for n, v in d.items()} for k, d in per_mg.items()}, '| break-even uses', round(break_even, 2),
      '| test GWP', round(test_gwp, 4))
