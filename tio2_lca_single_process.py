"""
Cradle-to-Gate LCA -- Single-Process Assessment (full-spectrum LED-assisted photothermal synthesis)
===================================================================================================
Baseline, non-comparative, cradle-to-gate LCA of the developed process only: functional unit (FU) = 1 g of
Cu/carbon-modified TiO2 catalyst produced. Foreground inventory as in Methods Section 2.15.2 (glucose 1.570 g,
CuSO4.5H2O 0.072 g, DI water 2.000 L, wastewater 2.000 L, centrifugation 10 min, drying 70 degC 12 h; electricity itemized per unit
operation, 10.616 kWh per FU). No comparator route is modeled.
This script computes the baseline only; the one-at-a-time sensitivity analysis (drying-oven electricity +/-20%,
post-wash wastewater 150/200 mL, catalyst yield 95%; manuscript Table 9) is computed by
tio2_lca_sensitivity_revision.py, which imports this module. No Monte Carlo propagation is performed.
Characterization factors: project-specific Ecoinvent 3 (cut-off) values embedded below; transport and packaging
use literature-average proxy factors for global-warming potential only, flagged as approximate (Section
"TRANSPORT AND PACKAGING" below).

Run: python3 tio2_lca_single_process.py
Writes: ../data/lca_single_process_results.json
"""

import json

# ---------------------------------------------------------------------------
# 1. CHARACTERIZATION FACTORS (author-supplied, Ecoinvent 3 - allocation,
#    cut-off), embedded directly. Values are per unit of the flow named
#    (kWh for electricity, kg for solids/water, m3 for wastewater).
# ---------------------------------------------------------------------------
MIDPOINT_CODES = ["GWP","SOD","IR","OF-HH","FPMF","OF-TE","TA","FE","ME",
                  "TET","FET","MET","HCT","HNT","LU","MRS","FRS","WC"]
MIDPOINT_UNITS = {
    "GWP": "kg CO2 eq", "SOD": "kg CFC-11 eq", "IR": "kBq Co-60 eq",
    "OF-HH": "kg NOx eq", "FPMF": "kg PM2.5 eq", "OF-TE": "kg NOx eq",
    "TA": "kg SO2 eq", "FE": "kg P eq", "ME": "kg N eq",
    "TET": "kg 1,4-DCB", "FET": "kg 1,4-DCB", "MET": "kg 1,4-DCB",
    "HCT": "kg 1,4-DCB", "HNT": "kg 1,4-DCB", "LU": "m2 x yr crop eq",
    "MRS": "kg Cu eq", "FRS": "kg oil eq", "WC": "m3",
}
MIDPOINT_NAMES = {
    "GWP": "Global warming potential", "SOD": "Stratospheric ozone depletion",
    "IR": "Ionizing radiation", "OF-HH": "Ozone formation, human health",
    "FPMF": "Fine particulate matter formation",
    "OF-TE": "Ozone formation, terrestrial ecosystems",
    "TA": "Terrestrial acidification", "FE": "Freshwater eutrophication",
    "ME": "Marine eutrophication", "TET": "Terrestrial ecotoxicity",
    "FET": "Freshwater ecotoxicity", "MET": "Marine ecotoxicity",
    "HCT": "Human carcinogenic toxicity", "HNT": "Human non-carcinogenic toxicity",
    "LU": "Land use", "MRS": "Mineral resource scarcity",
    "FRS": "Fossil resource scarcity", "WC": "Water consumption",
}

CF = {
  "Electricity_grid_TH": {
    "GWP": 0.797, "SOD": 2.13e-07, "IR": 0.00245, "OF-HH": 0.00127,
    "FPMF": 0.000732, "OF-TE": 0.00133, "TA": 0.0019, "FE": 0.00051,
    "ME": 3.25e-05, "TET": 2.1, "FET": 0.0492, "MET": 0.063, "HCT": 0.0383,
    "HNT": 0.851, "LU": 0.00489, "MRS": 0.00139, "FRS": 0.231, "WC": 0.00337,
    "DALY": 1.53e-06, "species_yr": 3.29e-09, "USD2013": 0.0549,
    "mPt_HH": 25.5, "mPt_ES": 0.889, "mPt_RS": 0.391,
  },
  "TiO2": {
    "GWP": 6.84, "SOD": 1.98e-06, "IR": 0.215, "OF-HH": 0.0179,
    "FPMF": 0.0159, "OF-TE": 0.0184, "TA": 0.041, "FE": 0.00427,
    "ME": 0.000266, "TET": 37.0, "FET": 0.482, "MET": 0.64, "HCT": 0.339,
    "HNT": 10.4, "LU": 0.321, "MRS": 0.97, "FRS": 1.6, "WC": 0.0659,
    "DALY": 2e-05, "species_yr": 3.75e-08, "USD2013": 0.603,
    "mPt_HH": 333.0, "mPt_ES": 10.1, "mPt_RS": 4.3,
  },
  "Glucose": {
    "GWP": 1.73, "SOD": 7.85e-06, "IR": 0.0827, "OF-HH": 0.00502,
    "FPMF": 0.00305, "OF-TE": 0.0052, "TA": 0.00873, "FE": 0.000661,
    "ME": 0.000976, "TET": 10.6, "FET": 0.106, "MET": 0.133, "HCT": 0.112,
    "HNT": 1.99, "LU": 0.478, "MRS": 0.0092, "FRS": 0.403, "WC": 0.122,
    "DALY": 4.59e-06, "species_yr": 1.37e-08, "USD2013": 0.0973,
    "mPt_HH": 76.4, "mPt_ES": 3.71, "mPt_RS": 0.694,
  },
  "CuSO4_5H2O": {
    "GWP": 1.99, "SOD": 3.92e-06, "IR": 0.112, "OF-HH": 0.018,
    "FPMF": 0.0385, "OF-TE": 0.0184, "TA": 0.108, "FE": 0.0169,
    "ME": 0.000383, "TET": 777.0, "FET": 14.2, "MET": 17.9, "HCT": 0.832,
    "HNT": 174.0, "LU": 0.13, "MRS": 0.285, "FRS": 0.57, "WC": 0.0416,
    "DALY": 6.86e-05, "species_yr": 6.45e-08, "USD2013": 0.201,
    "mPt_HH": 1140.0, "mPt_ES": 17.4, "mPt_RS": 1.43,
  },
  "DI_water": {
    "GWP": 0.000484, "SOD": 4.38e-10, "IR": 3.21e-05, "OF-HH": 1.14e-06,
    "FPMF": 1.25e-06, "OF-TE": 1.18e-06, "TA": 3.13e-06, "FE": 2.17e-07,
    "ME": 1.69e-08, "TET": 0.00717, "FET": 7.05e-05, "MET": 9.23e-05,
    "HCT": 5.41e-05, "HNT": 0.0014, "LU": 1.2e-05, "MRS": 5.81e-06,
    "FRS": 0.000114, "WC": 0.000965,
    "DALY": 3.58e-09, "species_yr": 1.43e-11, "USD2013": 2.83e-05,
    "mPt_HH": 0.0596, "mPt_ES": 0.00387, "mPt_RS": 0.000202,
  },
  "Wastewater": {
    "GWP": 0.286, "SOD": 8.23e-07, "IR": 0.0106, "OF-HH": 0.000682,
    "FPMF": 0.000426, "OF-TE": 0.000703, "TA": 0.000832, "FE": 0.00635,
    "ME": 0.00226, "TET": 1.21, "FET": 0.0391, "MET": 0.0531, "HCT": 0.0492,
    "HNT": 2.2, "LU": 0.00575, "MRS": 0.00219, "FRS": 0.0579, "WC": -1.9,
    "DALY": -2.43061e-06, "species_yr": -1.77184e-08, "USD2013": 0.014179226,
    "mPt_HH": -40.50679684, "mPt_ES": -4.791316609, "mPt_RS": 0.101159839,
  },
}

# ---------------------------------------------------------------------------
# 2. FOREGROUND INVENTORY, functional unit (FU) = 1 g of catalyst produced
#    (production/synthesis step only; MB-degradation, performance testing, and
#    wastewater-treatment use are out of scope). Quantities follow the Methods
#    of the manuscript and the validated recipe (author decision 2026-09-24/25):
#    1.570 g glucose, 0.072 g CuSO4.5H2O, 0.1 L + 2.0 L DI water, centrifugation
#    10 min, drying 70 degC for 12 h at 800 W.
# ---------------------------------------------------------------------------
ENERGY_BREAKDOWN_KWH = {
    "Magnetic stirrer (suspension preparation, 300 W, 30 min)": 0.150,
    "Ultrasonication (homogenization, 600 W, 30 min)": 0.300,
    "Full-spectrum LED (photothermal synthesis, 100 W, 140 min)": 0.233,
    "Centrifugation (catalyst separation, 2000 W, 10 min)": 2.0*10/60,
    "Drying oven (70 degC, 800 W, 12 h)": 9.600,
}
ELECTRICITY_TOTAL_KWH = sum(ENERGY_BREAKDOWN_KWH.values())  # 10.616 kWh / g

FOREGROUND = {
    "TiO2":               {"qty": 1.0,      "unit": "g"},
    "Glucose":            {"qty": 1.57,    "unit": "g"},
    "CuSO4_5H2O":         {"qty": 0.072,    "unit": "g"},
    "DI_water":           {"qty": 2.0,      "unit": "kg"},   # 0.100 L synthesis + 1.900 L washing, rho = 1 kg/L
    "Electricity_grid_TH":{"qty": ELECTRICITY_TOTAL_KWH, "unit": "kWh"},
    "Wastewater":         {"qty": 0.002000, "unit": "m3"},   # 50 mL process liquor + 1.950 L post-wash water (1.900 L wash + 50 mL synthesis water retained in the wet solid) = DI input
}
# unit conversion to the CF's native basis (kg for solids/water-mass flows; kWh; m3)
FOREGROUND_KG = {"TiO2": 0.001, "Glucose": 0.00157, "CuSO4_5H2O": 0.000072}  # g -> kg

# ---------------------------------------------------------------------------
# 3. TRANSPORT AND PACKAGING -- approximate/proxy characterization, GWP only
#    (author-supplied masses/distances; no project-specific Ecoinvent CF was
#    available for a generic "transport, freight, lorry" or "packaging,
#    plastic" process, so literature-average GWP factors are used as an
#    order-of-magnitude proxy; NOT extended to the other 17 midpoint or the
#    endpoint categories, both because no defensible proxy value was at hand
#    for those categories and because -- given the tiny masses/distances
#    involved -- their contribution would in any case be negligible relative
#    to electricity, TiO2, glucose, and CuSO4*5H2O, which is demonstrated
#    quantitatively for GWP below).
# ---------------------------------------------------------------------------
TRANSPORT_GWP_FACTOR_KG_CO2E_PER_TKM = 0.12   # generic lorry freight, literature-average proxy
PACKAGING_GWP_FACTOR_KG_CO2E_PER_KG = 2.5     # generic plastic (PP/PE) packaging production, literature-average proxy

raw_material_transport_tkm = 0.002 / 1000.0 * 500   # 0.002 kg over 500 km
product_transport_tkm      = 0.001 / 1000.0 * 100   # 0.001 kg over 100 km
packaging_kg = 0.005

transport_gwp = (raw_material_transport_tkm + product_transport_tkm) * TRANSPORT_GWP_FACTOR_KG_CO2E_PER_TKM
packaging_gwp = packaging_kg * PACKAGING_GWP_FACTOR_KG_CO2E_PER_KG
transport_packaging_gwp = transport_gwp + packaging_gwp

# ---------------------------------------------------------------------------
# 4. IMPACT ASSESSMENT
# ---------------------------------------------------------------------------
def flow_amount_in_cf_units(flow):
    if flow in FOREGROUND_KG:
        return FOREGROUND_KG[flow]
    return FOREGROUND[flow]["qty"]

midpoint_by_flow = {}
midpoint_total = {c: 0.0 for c in MIDPOINT_CODES}
for flow in FOREGROUND:
    amt = flow_amount_in_cf_units(flow)
    row = {}
    for code in MIDPOINT_CODES:
        v = amt * CF[flow][code]
        row[code] = v
        midpoint_total[code] += v
    midpoint_by_flow[flow] = row

# add transport+packaging GWP-only contribution
midpoint_total["GWP"] += transport_packaging_gwp
midpoint_by_flow["Transport_and_packaging"] = {c: (transport_packaging_gwp if c == "GWP" else None) for c in MIDPOINT_CODES}

endpoint_total = {"DALY": 0.0, "species_yr": 0.0, "USD2013": 0.0}
singlescore_total = {"mPt_HH": 0.0, "mPt_ES": 0.0, "mPt_RS": 0.0}
endpoint_by_flow = {}
singlescore_by_flow = {}
for flow in FOREGROUND:
    amt = flow_amount_in_cf_units(flow)
    e_row = {k: amt * CF[flow][k] for k in endpoint_total}
    s_row = {k: amt * CF[flow][k] for k in singlescore_total}
    endpoint_by_flow[flow] = e_row
    singlescore_by_flow[flow] = s_row
    for k in endpoint_total:
        endpoint_total[k] += e_row[k]
    for k in singlescore_total:
        singlescore_total[k] += s_row[k]
singlescore_total_mPt = sum(singlescore_total.values())

# GWP contribution shares (%), including transport+packaging, for the
# narrative/limitations discussion of how negligible that proxy addition is
gwp_shares_pct = {}
gwp_grand_total = midpoint_total["GWP"]
for flow in FOREGROUND:
    gwp_shares_pct[flow] = 100.0 * midpoint_by_flow[flow]["GWP"] / gwp_grand_total
gwp_shares_pct["Transport_and_packaging"] = 100.0 * transport_packaging_gwp / gwp_grand_total

# energy-equipment shares (%), for the process-hotspot discussion
energy_shares_pct = {k: 100.0 * v / ELECTRICITY_TOTAL_KWH for k, v in ENERGY_BREAKDOWN_KWH.items()}

# ---------------------------------------------------------------------------
# 5. ASSEMBLE OUTPUT
# ---------------------------------------------------------------------------
out = {
    "functional_unit": "1 g of Cu/carbon-modified TiO2 catalyst produced (full-spectrum LED-assisted photothermal synthesis; production/synthesis stage only)",
    "foreground_inventory": {
        "TiO2_g": 1.0, "Glucose_g": 1.57, "CuSO4_5H2O_g": 0.072,
        "DI_water_kg": 2.0, "Electricity_kWh": ELECTRICITY_TOTAL_KWH,
        "Wastewater_m3": 0.002000,
        "energy_breakdown_kWh": ENERGY_BREAKDOWN_KWH,
        "energy_shares_pct": energy_shares_pct,
    },
    "transport_and_packaging": {
        "raw_material_transport_tkm": raw_material_transport_tkm,
        "product_transport_tkm": product_transport_tkm,
        "packaging_kg": packaging_kg,
        "gwp_factor_transport_kgCO2e_per_tkm": TRANSPORT_GWP_FACTOR_KG_CO2E_PER_TKM,
        "gwp_factor_packaging_kgCO2e_per_kg": PACKAGING_GWP_FACTOR_KG_CO2E_PER_KG,
        "gwp_contribution_kgCO2e": transport_packaging_gwp,
        "gwp_contribution_pct_of_total": gwp_shares_pct["Transport_and_packaging"],
        "note": "approximate/proxy literature-average factors; GWP only; not extended to other categories (see script docstring)",
    },
    "midpoint_codes": MIDPOINT_CODES,
    "midpoint_units": MIDPOINT_UNITS,
    "midpoint_names": MIDPOINT_NAMES,
    "midpoint_total_per_FU": midpoint_total,
    "midpoint_by_flow": midpoint_by_flow,
    "gwp_contribution_shares_pct": gwp_shares_pct,
    "endpoint_total_per_FU": endpoint_total,
    "endpoint_by_flow": endpoint_by_flow,
    "singlescore_total_per_FU_mPt": singlescore_total,
    "singlescore_total_mPt_all": singlescore_total_mPt,
    "singlescore_by_flow": singlescore_by_flow,
}

with open("../data/lca_single_process_results.json", "w", encoding="utf-8") as f:
    json.dump(out, f, indent=2)

print("GWP total (kg CO2-eq / g catalyst):", round(midpoint_total["GWP"], 4))
print("  of which transport+packaging (proxy, GWP only): %.2f%%" % gwp_shares_pct["Transport_and_packaging"])
print("Water consumption (m3/g catalyst):", midpoint_total["WC"])
print("Single score total (mPt/g catalyst):", round(singlescore_total_mPt, 2))
print("Single score breakdown: HH=%.2f ES=%.2f RS=%.2f mPt" % (
    singlescore_total["mPt_HH"], singlescore_total["mPt_ES"], singlescore_total["mPt_RS"]))
print()
print("Energy breakdown (kWh / g catalyst, %% of total):")
for k, v in ENERGY_BREAKDOWN_KWH.items():
    print(f"  {k}: {v:.3f} kWh ({energy_shares_pct[k]:.1f}%)")
print()
print("GWP contribution by foreground flow (%):")
for k, v in gwp_shares_pct.items():
    print(f"  {k}: {v:.2f}%")
print("DONE -- wrote lca_single_process_results.json")
