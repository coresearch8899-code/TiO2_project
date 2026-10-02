# Reproducibility package: Cu/carbon-modified TiO₂, full-spectrum LED-assisted photothermal synthesis

This package regenerates every number and figure reported in the manuscript and the Supplementary Material.

## Run

```
cd code
python run_all.py
```

`run_all.py` runs the 15 scripts below in dependency order and stops at the first failure. A full run takes about 2 minutes; the ML benchmark and the figures take most of that time.

Requirements: Python 3.11 and the package versions pinned in `requirements.txt` (`pip install -r requirements.txt`). These are the versions of the verified run; other versions of xgboost, scikit-learn or scipy can change the last significant figures of the LOOCV and SLSQP values.

## Layout

| Folder | Contents |
|---|---|
| `code/` | analysis and figure scripts (this folder) |
| `data/` | raw laboratory data (absorbance workbook `data_analyse_1.xlsx`, diffuse reflectance `%R.xlsx`, `RAW/`, `SEM TiO2/`, `FTIR TiO2 N7-N12/`, `XRD-26.08.13/`, …) and the JSON result files written by the scripts |
| `figures/` | figure output, 600 dpi PNG plus vector PDF |

`data/TiO2-SolarModified_LCA_v4.xlsx` is the original laboratory source workbook for the LCA characterization factors embedded in `tio2_lca_single_process.py`. Its filename is kept unchanged for provenance. The factors for Thailand-grid electricity, TiO₂, D-glucose, CuSO₄·5H₂O, deionized water and wastewater treatment come from it. Its rooftop-photovoltaic row is not used, and no photovoltaic or comparative solar-energy scenario is modeled.

## Scripts, outputs and where they are reported

| Step | Script | Writes | Reported in |
|---|---|---|---|
| 1 | `tio2_reanalysis.py` | `data/manuscript_data.json`, `code/reanalysis_outputs_v2/*.csv` | Tables 2, 4, 6, 7; Tables S1–S3, S6, S7, S9–S14 |
| 2 | `tio2_reoptimize_hull.py` | `data/optimization_hull_results.json`; updates `manuscript_data.json` | Table 5 (including the D = 1 region), Table S15 |
| 3 | `tio2_synergy_analysis.py` | `data/synergy_analysis_results.json` | Section 3.6, Table S5 |
| 4 | `tio2_dark_adsorption_analysis.py` | `data/dark_adsorption_results.json` | Table 1, Tables S8a–c |
| 5 | `tio2_lca_single_process.py` | `data/lca_single_process_results.json` | Table 8, Tables S16–S18 |
| 6 | `tio2_lca_sensitivity_revision.py` | `data/lca_sensitivity_revision_results.json` | Table 9 |
| 7 | `tio2_antibacterial_analysis.py` | `data/antibacterial_results.json` | Table S20 |
| 8 | `tio2_characterization_figures.py` | `fig01_…`, `fig02_…` | Figs. 1–2 |
| 9 | `tio2_figures.py` | `fig03_…`, `fig04_…`, `fig05_…` | Figs. 3–5 |
| 10–13 | `fig06_ml_benchmark.py`, `fig07_pareto_hull.py`, `fig08_validation.py`, `fig09_ablation_and_graphical_abstract.py` | `fig06_…` to `fig09_…`, `fig00_graphical_abstract` | Figs. 6–9, graphical abstract |
| 14 | `tio2_lca_single_process_figure.py` | `fig10_lca_single_process` | Fig. 10 |
| 15 | `tio2_antibacterial_figures.py` | `figS1_antibacterial_kinetics` | Fig. S1 |

`tio2_reoptimize_hull.py --check` first reproduces the earlier spherical-bound optimization (‖x‖₂ ≤ √2) and stops.

## Numerical reproducibility

A run on 2026-09-25 reproduced every figure pixel for pixel. It also reproduced every JSON value to the precision reported in the manuscript (3 decimals).

Two groups of values differ from run to run in the last few significant figures, well below the reported precision:

- the XGBoost leave-one-out values, around the 8th significant figure, due to floating-point summation order in multithreaded XGBoost;
- the SLSQP bounds of the D = 1 region, below 0.001 in every reported unit.
