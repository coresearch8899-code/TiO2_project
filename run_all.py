"""Regenerate every reported number and figure, in dependency order.

Run from this folder:  python run_all.py
Inputs: ../data (raw laboratory data; raw absorbance and plate counts are also embedded in the scripts).
Outputs: ../data/*.json, ./reanalysis_outputs_v2/*.csv, ../figures/*.png|pdf. See README.md for the file-to-table map.
"""
import subprocess, sys, time

STEPS = [
    # analyses
    ('tio2_reanalysis.py', 'RSM, ML benchmark, spherical-bound optimization, ablation, validation -> ../data/manuscript_data.json'),
    ('tio2_reoptimize_hull.py', 'convex-hull re-optimization -> Table 5, Table S15; updates manuscript_data.json'),
    ('tio2_synergy_analysis.py', 'Cu x carbon interaction test -> Section 3.6, Table S5'),
    ('tio2_dark_adsorption_analysis.py', 'dark adsorption / total removal -> Table 1, Tables S8a-c'),
    ('tio2_lca_single_process.py', 'cradle-to-gate LCA baseline -> Table 8, Tables S16-S18'),
    ('tio2_lca_sensitivity_revision.py', 'one-at-a-time LCA sensitivity -> Table 9'),
    ('tio2_antibacterial_analysis.py', 'antibacterial kinetics -> Table S20'),
    # figures
    ('tio2_characterization_figures.py', 'Figs. 1-2'),
    ('tio2_figures.py', 'Figs. 3-5'),
    ('fig06_ml_benchmark.py', 'Fig. 6'),
    ('fig07_pareto_hull.py', 'Fig. 7'),
    ('fig08_validation.py', 'Fig. 8'),
    ('fig09_ablation_and_graphical_abstract.py', 'Fig. 9 and graphical abstract'),
    ('tio2_lca_single_process_figure.py', 'Fig. 10'),
    ('tio2_antibacterial_figures.py', 'Fig. S1'),
]

for script, what in STEPS:
    t0 = time.time()
    print(f'\n=== {script}  ({what})', flush=True)
    r = subprocess.run([sys.executable, script], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding='utf-8', errors='replace')
    if r.returncode:
        print(r.stdout[-4000:])
        sys.exit(f'FAILED: {script}')
    print(f'    ok ({time.time() - t0:.0f} s)')
print('\nAll steps completed.')
