"""Regenerate every reported number and figure, in the order of the analysis.

Run from this folder:  python run_all.py
Inputs: ../data (raw laboratory data; raw absorbance and plate counts are also embedded in the scripts).
Outputs: ../data/*.json, ./reanalysis_outputs_v2/*.csv, ../figures/*.png|pdf. See ../README.md for the file-to-table map.
"""
import subprocess, sys, time

STEPS = [
    # characterization of the surface-modified TiO2
    (['s01_characterization_figures.py'], 'XRD, FTIR, SEM, DRS/Tauc figures'),
    # process definition: response surfaces, ANOVA, ML benchmark, validation, ablation kinetics
    (['s02_rsm_ml_validation_kinetics.py'], 'RSM, ANOVA, LOOCV, ML benchmark, validation, ablation kinetics -> ../data/manuscript_data.json'),
    (['s03_lack_of_fit_design_hull.py'], 'lack of fit; convex hull of the tested points -> ../data/design_support_results.json'),
    (['s04_dark_adsorption.py'], 'dark adsorption / total removal'),
    (['s05_rsm_surface_figures.py'], 'response-surface figures'),
    (['s06_ml_benchmark_figure.py'], 'ML benchmark figure'),
    (['s07_optimization_tested_hull.py', '--bbd'], 'optimization over the nominal BBD region -> ../data/optimization_hull_results.json'),
    (['s07_optimization_tested_hull.py'], 'optimization within the tested hull -> ../data/optimization_tested_hull_results.json'),
    (['s08_pareto_figure.py'], 'Pareto figure'),
    (['s09_validation_figure.py'], 'validation figure'),
    (['s10_synergy_test.py'], 'Cu x carbon interaction test'),
    # life cycle assessment of the validated route
    (['s11_lca_baseline.py'], 'cradle-to-gate LCA baseline -> ../data/lca_single_process_results.json'),
    (['s12_lca_sensitivity.py'], 'one-at-a-time sensitivity -> ../data/lca_sensitivity_revision_results.json'),
    (['s13_lca_scenarios.py'], 'process-improvement scenarios -> ../data/lca_scenario_results.json'),
    (['s14_ablation_figure_graphical_abstract.py'], 'ablation figure and graphical abstract (needs s04 and s11 results)'),
    (['s15_lca_hotspot_figure.py'], 'LCA hotspot figure'),
    (['s16_lca_scenario_figure.py'], 'LCA scenario figure'),
    (['s17_system_boundary_figure.py'], 'system boundary figure'),
    # supplementary antibacterial characterization
    (['s18_antibacterial_kinetics.py'], 'antibacterial kinetics'),
    (['s19_antibacterial_figure.py'], 'antibacterial figure'),
    (['s20_graphical_abstract.py'], 'graphical abstract (ET&I)'),
]

for args, what in STEPS:
    t0 = time.time()
    print(f"\n=== {' '.join(args)}  ({what})", flush=True)
    r = subprocess.run([sys.executable, *args], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding='utf-8', errors='replace')
    if r.returncode:
        print(r.stdout[-4000:])
        sys.exit(f"FAILED: {' '.join(args)}")
    print(f'    ok ({time.time() - t0:.0f} s)')
print('\nAll steps completed.')
