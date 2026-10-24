"""Regenerate every reported number and figure, in the order of the analysis (manuscript Sections 3.1-3.8).

Run from this folder:  python run_all.py
Inputs: ../data (raw laboratory data; raw absorbance and plate counts are also embedded in the scripts).
Outputs: ../data/*.json, ./reanalysis_outputs_v2/*.csv, ../figures/*.png|pdf. See ../README.md for the file-to-table map.
"""
import subprocess, sys, time

STEPS = [
    # 3.1 physicochemical characterization
    (['s01_characterization_figures.py'], 'XRD, FTIR, SEM, DRS/Tauc -> Figs. 1-2'),
    # 3.2-3.3 response surfaces, ANOVA, ML benchmark; validation and ablation kinetics are computed here too
    (['s02_rsm_ml_validation_kinetics.py'], 'RSM, ANOVA, LOOCV, ML benchmark, validation, ablation kinetics -> ../data/manuscript_data.json'),
    (['s03_lack_of_fit_design_hull.py'], 'lack of fit; convex hull of the tested points -> ../data/design_support_results.json'),
    (['s04_dark_adsorption.py'], 'dark adsorption / total removal -> Table 1, Table S8'),
    (['s05_fig03_05_rsm_surfaces.py'], 'Figs. 3-5'),
    (['s06_fig06_ml_benchmark.py'], 'Fig. 6'),
    # 3.4 multi-objective optimization (full-BBD region first, kept for comparison; then the tested hull)
    (['s07_optimization_tested_hull.py', '--bbd'], 'optimization over the nominal BBD region -> ../data/optimization_hull_results.json'),
    (['s07_optimization_tested_hull.py'], 'optimization within the tested hull -> Table 5, Table S15'),
    (['s08_fig07_pareto.py'], 'Fig. 7'),
    # 3.5 experimental validation
    (['s09_fig08_validation.py'], 'Fig. 8'),
    # 3.6 mechanistic ablation
    (['s10_synergy_test.py'], 'Cu x carbon interaction test -> Section 3.6, Table S5'),
    # 3.7 life cycle assessment
    (['s11_lca_baseline.py'], 'cradle-to-gate LCA -> Table 8, Tables S16-S18'),
    (['s12_lca_sensitivity.py'], 'one-at-a-time sensitivity -> Table 9'),
    (['s13_fig09_ablation_graphical_abstract.py'], 'Fig. 9 and graphical abstract (needs s04 and s11 results)'),
    (['s14_fig10_lca.py'], 'Fig. 10'),
    # Supplementary Material S8: antibacterial characterization
    (['s15_antibacterial_kinetics.py'], 'antibacterial kinetics -> Table S20'),
    (['s16_figS1_antibacterial.py'], 'Fig. S1'),
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
