"""
TiO2 Cu/Carbon-Modified Photocatalyst -- Full Re-Analysis
===================================================================================
Canonical, self-contained re-analysis script. Raw absorbance data as supplied in
data_analyse_1.xlsx (2026-09-20 update, third data revision):
  - 15-run Box-Behnken absorbance data unchanged from the previous (2026-09-19)
    correction.
  - Apparent quantum efficiency (AQE) computed from a directly measured
    photosynthetic photon flux density (PPFD = 1361.706 umol photons m^-2 s^-1)
    rather than derived from peak-wavelength photon energy (Planck's equation) --
    this removes the monochromatic-source approximation used in an earlier round.
  - Specific Energy Consumption (SEC, Wh per mg MB photodegraded) is carried as a
    5th formal response alongside Photodeg%, k_app, AQE, and EEO, split into a
    synthesis-energy share (SEC_syn) and a photocatalytic-test-energy share
    (SEC_photo), with SEC_total = SEC_syn + SEC_photo.
  - The constrained-optimization-identified operating condition supplied for
    experimental validation has been revised to Glucose=1.57 g, CuSO4=0.072 g,
    Time=140 min (previously 1.505 g / 0.072 g / 135.9 min). V1-V3 are direct
    triplicate replicates of this new point; V4-V7, previously dropped for lack
    of data, are now supplied as single-factor (and joint glucose+time)
    perturbations around it and are analyzed as a local sensitivity/robustness
    check, each compared against the RSM model's own prediction at that set's
    actual tested condition (not the fixed new-optimum prediction).
  - Mechanistic ablation kinetics (K1,K2,K4,K5,K6) re-run at the same new
    condition (1.57 g / 0.072 g / 140 min); K3 (bare TiO2, no synthesis
    modification at all: 0 g / 0 g / 0 min) is the untreated commercial-TiO2
    control, as supplied. K6 glucose = 1.57 g, confirmed by the laboratory
    (2026-09-25); the entry of 157 g in data_analyse_1.xlsx is a typographical error.

Sections
--------
  1. Constants and embedded raw data
  2. Response functions (dark-adsorption-free; AQE via measured PPFD; SEC added)
  3. Box-Behnken response table + RSM quadratic model fitting (5 responses)
  4. Machine-learning surrogate benchmarking (9 algorithms, in-sample + LOOCV, 5 responses)
  5. Multi-objective optimization (3 objectives: Photodeg% max, EEO min, SEC min)
  6. Mechanistic ablation kinetics (K1-K6)
  7. Experimental validation (V1-V7) vs. (a) the RSM model's own prediction at
     each set's actual tested condition, (b) the original design-center
     condition (Runs 13-15), and (c) the new-optimum replicate mean (V1-V3)
  8. Export ../data/manuscript_data.json (merging in the untouched characterization /
     antibacterial / convex-hull optimization sections of the existing file); the
     life cycle assessment is computed separately (tio2_lca_single_process.py)
"""

import os, json
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import differential_evolution
import statsmodels.api as sm
import statsmodels.formula.api as smf

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel as C
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import r2_score, mean_squared_error
import xgboost as xgb
import lightgbm as lgb

OUTDIR = "reanalysis_outputs_v2"
os.makedirs(OUTDIR, exist_ok=True)
RNG = 42

# ---------------------------------------------------------------------------
# 1. CONSTANTS AND EMBEDDED RAW DATA (source: data_analyse_1.xlsx)
# ---------------------------------------------------------------------------

LED_POWER_W = 100.0
BEAKER_AREA_M2 = 0.00196
POWER_DENSITY_W_M2 = 282.857
PEAK_WAVELENGTH_NM = 619.4
PPFD_UMOL_M2_S = 1361.706            # NEW: directly measured photon flux density
MB_MOLAR_MASS_G_MOL = 319.85
SAMPLE_VOLUME_L = 0.1
INITIAL_MB_CONC_MG_L = 2.5
INITIAL_MB_MASS_MG = INITIAL_MB_CONC_MG_L * SAMPLE_VOLUME_L   # 0.25 mg
PHOTO_TIME_MIN = 60.0

H_PLANCK = 6.62607015e-34
C_LIGHT = 2.99792458e8
N_AVOGADRO = 6.02214076e23

def moles_incident_photons_ppfd(time_min=PHOTO_TIME_MIN):
    """Photon count from directly measured PPFD (umol photons m^-2 s^-1)."""
    mol_photons_per_s = PPFD_UMOL_M2_S * 1e-6 * BEAKER_AREA_M2
    return mol_photons_per_s * (time_min * 60.0)

N_PHOTONS_MOL_60MIN = moles_incident_photons_ppfd(60.0)

# --- Box-Behnken design matrix, corrected (Run, Glucose_g, CuSO4_g, Time_min, Ai, A0, A60) ---
BBD_RAW = pd.DataFrame([
    [1,  0.5,  0.02, 150, 0.571, 0.498, 0.419],
    [2,  2.0,  0.02, 150, 0.627, 0.514, 0.492],
    [3,  0.5,  0.08, 150, 0.571, 0.474, 0.447],
    [4,  2.0,  0.05, 150, 0.571, 0.475, 0.402],
    [5,  0.5,  0.05, 120, 0.571, 0.464, 0.398],
    [6,  2.0,  0.05, 120, 0.571, 0.529, 0.457],
    [7,  0.5,  0.05, 180, 0.571, 0.463, 0.421],
    [8,  2.0,  0.05, 180, 0.571, 0.499, 0.416],
    [9,  1.25, 0.02, 120, 0.627, 0.527, 0.451],
    [10, 1.25, 0.08, 120, 0.627, 0.550, 0.434],
    [11, 1.25, 0.02, 180, 0.627, 0.541, 0.501],
    [12, 1.25, 0.08, 180, 0.627, 0.551, 0.498],
    [13, 1.25, 0.05, 150, 0.627, 0.511, 0.397],
    [14, 1.25, 0.05, 150, 0.627, 0.513, 0.399],
    [15, 1.25, 0.05, 150, 0.627, 0.514, 0.401],
], columns=["Run", "Glucose_g", "CuSO4_g", "Time_min", "Ai", "A0", "A60"])

# --- Validation, new optimum condition (Glucose=1.57 g, CuSO4=0.072 g, Time=140 min) ---
# V1-V3: direct replicates of the new optimum. V4-V5: joint glucose+time perturbation
# (-0.10 g / -5 min and +0.10 g / +5 min). V6-V7: CuSO4 perturbation (-0.012 g / +0.008 g).
VALIDATION_RAW = pd.DataFrame([
    ["V1", 1.57, 0.072, 140, 1, 0.587, 0.481,  0.371],
    ["V1", 1.57, 0.072, 140, 2, 0.587, 0.480,  0.370],
    ["V1", 1.57, 0.072, 140, 3, 0.587, 0.479,  0.369],
    ["V2", 1.57, 0.072, 140, 1, 0.587, 0.4812, 0.368],
    ["V2", 1.57, 0.072, 140, 2, 0.587, 0.4799, 0.375],
    ["V2", 1.57, 0.072, 140, 3, 0.587, 0.483,  0.378],
    ["V3", 1.57, 0.072, 140, 1, 0.587, 0.481,  0.376],
    ["V3", 1.57, 0.072, 140, 2, 0.587, 0.482,  0.380],
    ["V3", 1.57, 0.072, 140, 3, 0.587, 0.480,  0.375],
    ["V4", 1.47, 0.072, 135, 1, 0.587, 0.506,  0.432],
    ["V4", 1.47, 0.072, 135, 2, 0.587, 0.505,  0.432],
    ["V4", 1.47, 0.072, 135, 3, 0.587, 0.505,  0.432],
    ["V5", 1.67, 0.072, 145, 1, 0.587, 0.476,  0.406],
    ["V5", 1.67, 0.072, 145, 2, 0.587, 0.475,  0.406],
    ["V5", 1.67, 0.072, 145, 3, 0.587, 0.477,  0.406],
    ["V6", 1.57, 0.060, 140, 1, 0.587, 0.453,  0.386],
    ["V6", 1.57, 0.060, 140, 2, 0.587, 0.455,  0.387],
    ["V6", 1.57, 0.060, 140, 3, 0.587, 0.455,  0.386],
    ["V7", 1.57, 0.080, 140, 1, 0.587, 0.437,  0.370],
    ["V7", 1.57, 0.080, 140, 2, 0.587, 0.439,  0.371],
    ["V7", 1.57, 0.080, 140, 3, 0.587, 0.436,  0.372],
], columns=["Set", "Glucose_g", "CuSO4_g", "Time_min", "Rep", "Ai", "A0", "A60"])

VALIDATION_SET_LABELS = {
    "V1": "Direct replicate of new optimum", "V2": "Direct replicate of new optimum",
    "V3": "Direct replicate of new optimum",
    "V4": "Glucose -0.10 g, time -5 min", "V5": "Glucose +0.10 g, time +5 min",
    "V6": "CuSO4 -0.012 g", "V7": "CuSO4 +0.008 g",
}

# --- Mechanistic ablation kinetics K1,K2,K4,K5,K6 at the new condition; K3 = bare TiO2 (no synthesis modification) ---
KINETICS_RAW = pd.DataFrame([
    ["K1", 1.57, 0.072, 140, 1, 0.542, 0.466, 0.453, 0.427, 0.407, 0.379, "Heat / no light"],
    ["K1", 1.57, 0.072, 140, 2, 0.542333, 0.467, 0.450, 0.425, 0.408, 0.378, "Heat / no light"],
    ["K1", 1.57, 0.072, 140, 3, 0.542333, 0.465, 0.450, 0.424, 0.409, 0.378, "Heat / no light"],
    ["K2", 1.57, 0.072, 140, 1, 0.542333, 0.491, 0.460, 0.432, 0.419, 0.400, "Light / no added heat"],
    ["K2", 1.57, 0.072, 140, 2, 0.542333, 0.490, 0.459, 0.432, 0.418, 0.400, "Light / no added heat"],
    ["K2", 1.57, 0.072, 140, 3, 0.542333, 0.489, 0.459, 0.433, 0.417, 0.400, "Light / no added heat"],
    ["K3", 0.00, 0.00, 0,     1, 0.542333, 0.473, 0.451, 0.438, 0.416, 0.407, "Bare TiO2, light+heat"],
    ["K3", 0.00, 0.00, 0,     2, 0.542333, 0.474, 0.450, 0.438, 0.417, 0.409, "Bare TiO2, light+heat"],
    ["K3", 0.00, 0.00, 0,     3, 0.542333, 0.475, 0.459, 0.438, 0.417, 0.407, "Bare TiO2, light+heat"],
    ["K4", 1.57, 0.00, 140,   1, 0.542333, 0.461, 0.445, 0.421, 0.403, 0.377, "Carbon-only, light+heat"],
    ["K4", 1.57, 0.00, 140,   2, 0.542333, 0.462, 0.444, 0.420, 0.405, 0.370, "Carbon-only, light+heat"],
    ["K4", 1.57, 0.00, 140,   3, 0.542333, 0.462, 0.444, 0.420, 0.405, 0.375, "Carbon-only, light+heat"],
    ["K5", 0.00, 0.072, 140,  1, 0.542333, 0.477, 0.456, 0.435, 0.415, 0.384, "Cu-only, light+heat"],
    ["K5", 0.00, 0.072, 140,  2, 0.542333, 0.476, 0.456, 0.436, 0.415, 0.388, "Cu-only, light+heat"],
    ["K5", 0.00, 0.072, 140,  3, 0.542333, 0.478, 0.456, 0.436, 0.416, 0.388, "Cu-only, light+heat"],
    ["K6", 1.57, 0.072, 140,  1, 0.542333, 0.481, 0.464, 0.442, 0.407, 0.372, "Full Cu+C, light+heat"],
    ["K6", 1.57, 0.072, 140,  2, 0.542333, 0.482, 0.464, 0.441, 0.407, 0.372, "Full Cu+C, light+heat"],
    ["K6", 1.57, 0.072, 140,  3, 0.542333, 0.482, 0.464, 0.441, 0.408, 0.375, "Full Cu+C, light+heat"],
], columns=["Condition", "Glucose_g", "CuSO4_g", "Time_min", "Rep", "Ai", "A0",
            "A15", "A30", "A45", "A60", "Method"])


# ---------------------------------------------------------------------------
# 2. RESPONSE FUNCTIONS (dark-adsorption-free; SEC added)
# ---------------------------------------------------------------------------

def photodegradation_pct(a0, a60):
    return (a0 - a60) / a0 * 100.0

def k_app_per_min(a0, a60, t_min=PHOTO_TIME_MIN):
    return np.log(a0 / a60) / t_min

def mass_photodegraded_mg(a0, a60, initial_mass_mg=INITIAL_MB_MASS_MG):
    return initial_mass_mg * (a0 - a60) / a0

def aqe_pct(a0, a60, n_photons_mol=N_PHOTONS_MOL_60MIN):
    mass_mg = mass_photodegraded_mg(a0, a60)
    moles_mb = (mass_mg / 1000.0) / MB_MOLAR_MASS_G_MOL
    return moles_mb / n_photons_mol * 100.0

def eeo_kwh_m3_order(power_w, time_min, volume_l, a0, a60):
    p_kw = power_w / 1000.0
    t_h = time_min / 60.0
    log_ratio = np.log10(a0 / a60)
    return (p_kw * t_h * 1000.0) / (volume_l * log_ratio)

def synthesis_energy_wh(time_min, power_w=LED_POWER_W):
    return power_w * (time_min / 60.0)

def photo_energy_wh(power_w=LED_POWER_W, time_min=PHOTO_TIME_MIN):
    return power_w * (time_min / 60.0)

def sec_components_wh_per_mg(a0, a60, synth_time_min):
    """Returns (SEC_syn, SEC_photo, SEC_total), all Wh per mg MB photodegraded."""
    mass_mg = mass_photodegraded_mg(a0, a60)
    e_syn = synthesis_energy_wh(synth_time_min)
    e_photo = photo_energy_wh()
    sec_syn = e_syn / mass_mg
    sec_photo = e_photo / mass_mg
    sec_total = sec_syn + sec_photo
    return sec_syn, sec_photo, sec_total


def build_response_table(df):
    out = df.copy()
    out["Photodeg_pct"] = photodegradation_pct(out["A0"], out["A60"])
    out["k_app_min-1"] = k_app_per_min(out["A0"], out["A60"])
    out["Mass_photodeg_mg"] = mass_photodegraded_mg(out["A0"], out["A60"])
    out["AQE_pct"] = aqe_pct(out["A0"], out["A60"])
    out["EEO_kWh_m3_order"] = eeo_kwh_m3_order(LED_POWER_W, PHOTO_TIME_MIN,
                                                SAMPLE_VOLUME_L, out["A0"], out["A60"])
    sec_syn, sec_photo, sec_total = sec_components_wh_per_mg(out["A0"], out["A60"], out["Time_min"])
    out["SEC_syn_Wh_mg"] = sec_syn
    out["SEC_photo_Wh_mg"] = sec_photo
    out["SEC_total_Wh_mg"] = sec_total
    return out

print("=" * 78)
print("SECTION 1-2: RAW DATA (v2) LOADED, RESPONSES COMPUTED (incl. SEC, PPFD-based AQE)")
print("=" * 78)

bbd = build_response_table(BBD_RAW)
bbd.to_csv(f"{OUTDIR}/table1_bbd_v2.csv", index=False)
print(bbd[["Run","Glucose_g","CuSO4_g","Time_min","Photodeg_pct","k_app_min-1",
           "AQE_pct","EEO_kWh_m3_order","SEC_total_Wh_mg"]].round(4).to_string(index=False))

# ---------------------------------------------------------------------------
# 3. RSM MODELS + ANOVA (5 responses)
# ---------------------------------------------------------------------------

def code(x, lo, hi):
    mid = (lo + hi) / 2.0
    half = (hi - lo) / 2.0
    return (x - mid) / half

bbd["X1"] = code(bbd["Glucose_g"], 0.5, 2.0)
bbd["X2"] = code(bbd["CuSO4_g"], 0.02, 0.08)
bbd["X3"] = code(bbd["Time_min"], 120, 180)

RESPONSES = ["Photodeg_pct", "k_app_min-1", "AQE_pct", "EEO_kWh_m3_order", "SEC_total_Wh_mg"]
rsm_col_safe = {"k_app_min-1": "k_app", "Photodeg_pct": "Photodeg_pct",
                "AQE_pct": "AQE_pct", "EEO_kWh_m3_order": "EEO", "SEC_total_Wh_mg": "SEC_total"}

quad_formula = "{y} ~ X1 + X2 + X3 + I(X1**2) + I(X2**2) + I(X3**2) + X1:X2 + X1:X3 + X2:X3"

rsm_models = {}
rsm_summary_rows = []
anova_dict = {}
print("\n" + "=" * 78 + "\nSECTION 3: RSM MODELS\n" + "=" * 78)
for resp in RESPONSES:
    safe = rsm_col_safe[resp]
    d = bbd.rename(columns={resp: safe})
    model = smf.ols(quad_formula.format(y=safe), data=d).fit()
    rsm_models[resp] = model
    anova = sm.stats.anova_lm(model, typ=2)
    anova_dict[resp] = json.loads(anova.reset_index().rename(columns={"index": "Unnamed: 0"}).to_json(orient="records"))
    print(f"\n--- {resp} ---")
    print(f"R2 = {model.rsquared:.4f}   Adj-R2 = {model.rsquared_adj:.4f}   F-p = {model.f_pvalue:.4g}")
    rsm_summary_rows.append({
        "Response": resp, "R2": model.rsquared, "Adj_R2": model.rsquared_adj,
        "F_pvalue": model.f_pvalue, "RMSE": float(np.sqrt(model.mse_resid)),
    })

rsm_summary = pd.DataFrame(rsm_summary_rows)
print("\n", rsm_summary.round(4).to_string(index=False))

rsm_coefs = {resp: rsm_models[resp].params.to_dict() for resp in RESPONSES}

# ---------------------------------------------------------------------------
# 4. ML SURROGATE BENCHMARK (in-sample + LOOCV), 5 responses
# ---------------------------------------------------------------------------

print("\n" + "=" * 78 + "\nSECTION 4: ML BENCHMARK\n" + "=" * 78)
X_raw = bbd[["Glucose_g", "CuSO4_g", "Time_min"]].values

def make_models():
    return {
        "GPR": GaussianProcessRegressor(kernel=C(1.0) * RBF(length_scale=[1,1,1]) + WhiteKernel(1e-3),
                                         normalize_y=True, n_restarts_optimizer=5, random_state=RNG),
        "RF": RandomForestRegressor(n_estimators=300, random_state=RNG),
        "SVR": SVR(kernel="rbf", C=10.0, epsilon=0.01),
        "KNR": KNeighborsRegressor(n_neighbors=3),
        "GBR": GradientBoostingRegressor(random_state=RNG),
        "MLP": MLPRegressor(hidden_layer_sizes=(16,16), max_iter=5000, random_state=RNG),
        "XGBR": xgb.XGBRegressor(n_estimators=200, max_depth=3, random_state=RNG, verbosity=0),
        "LGBMR": lgb.LGBMRegressor(n_estimators=200, min_data_in_leaf=1, min_data_in_bin=1,
                                    verbosity=-1, random_state=RNG),
        "DTR": DecisionTreeRegressor(random_state=RNG),
    }

ml_rows = []
for resp in RESPONSES:
    y_raw = bbd[resp].values
    scaler_x, scaler_y = StandardScaler(), StandardScaler()
    Xs = scaler_x.fit_transform(X_raw)
    ys = scaler_y.fit_transform(y_raw.reshape(-1,1)).ravel()
    row = {"Response": resp}
    for name in make_models().keys():
        model = make_models()[name]
        model.fit(Xs, ys)
        r2_in = r2_score(ys, model.predict(Xs))
        loo = LeaveOneOut()
        y_true_cv, y_pred_cv = [], []
        for tr_idx, te_idx in loo.split(X_raw):
            sx, sy = StandardScaler(), StandardScaler()
            Xtr = sx.fit_transform(X_raw[tr_idx]); ytr = sy.fit_transform(y_raw[tr_idx].reshape(-1,1)).ravel()
            Xte = sx.transform(X_raw[te_idx])
            m = make_models()[name]; m.fit(Xtr, ytr)
            pred = sy.inverse_transform(m.predict(Xte).reshape(-1,1)).ravel()
            y_true_cv.append(y_raw[te_idx][0]); y_pred_cv.append(pred[0])
        r2_loo = r2_score(y_true_cv, y_pred_cv)
        rmse_loo = np.sqrt(mean_squared_error(y_true_cv, y_pred_cv))
        row[f"{name}_R2_insample"] = r2_in
        row[f"{name}_R2_LOOCV"] = r2_loo
        row[f"{name}_RMSE_LOOCV"] = rmse_loo
    # RSM LOOCV
    y_true_cv, y_pred_cv = [], []
    safe = rsm_col_safe[resp]
    for i in range(len(bbd)):
        train = bbd.drop(index=bbd.index[i]).rename(columns={resp: safe})
        test = bbd.iloc[[i]].rename(columns={resp: safe})
        m = smf.ols(quad_formula.format(y=safe), data=train).fit()
        pred = m.predict(test)[test.index[0]]
        y_true_cv.append(bbd.iloc[i][resp]); y_pred_cv.append(pred)
    row["RSM_R2_insample"] = rsm_models[resp].rsquared
    row["RSM_R2_LOOCV"] = r2_score(y_true_cv, y_pred_cv)
    row["RSM_RMSE_LOOCV"] = np.sqrt(mean_squared_error(y_true_cv, y_pred_cv))
    ml_rows.append(row)
    print(f"{resp}: RSM LOOCV R2={row['RSM_R2_LOOCV']:.3f}  (best ML LOOCV R2={max(row[f'{n}_R2_LOOCV'] for n in make_models()):.3f})")

ml_summary = pd.DataFrame(ml_rows)

# ---------------------------------------------------------------------------
# 5. MULTI-OBJECTIVE OPTIMIZATION (3 objectives: Photodeg% max, EEO min, SEC min)
# ---------------------------------------------------------------------------

print("\n" + "=" * 78 + "\nSECTION 5: MULTI-OBJECTIVE OPTIMIZATION (3 objectives)\n" + "=" * 78)

BOUNDS_REAL = [(0.5, 2.0), (0.02, 0.08), (120, 180)]
BBD_SUPPORT_RADIUS = np.sqrt(2.0)

_RSM_COEFS = {resp: rsm_models[resp].params.values.copy() for resp in RESPONSES}

def _poly_features(x1, x2, x3):
    ones = np.ones_like(np.asarray(x1, dtype=float))
    return np.column_stack([ones, x1, x2, x3, x1**2, x2**2, x3**2, x1*x2, x1*x3, x2*x3])

def rsm_predict(resp, glucose, cuso4, time):
    x1, x2, x3 = code(glucose,0.5,2.0), code(cuso4,0.02,0.08), code(time,120,180)
    feats = _poly_features(np.atleast_1d(x1), np.atleast_1d(x2), np.atleast_1d(x3))
    pred = feats @ _RSM_COEFS[resp]
    return float(pred[0]) if np.isscalar(glucose) or np.ndim(glucose)==0 else pred

def coded_norm(glucose, cuso4, time):
    x1, x2, x3 = code(glucose,0.5,2.0), code(cuso4,0.02,0.08), code(time,120,180)
    return np.sqrt(x1**2 + x2**2 + x3**2)

def support_penalty(glucose, cuso4, time, weight=50.0):
    excess = max(0.0, coded_norm(glucose, cuso4, time) - BBD_SUPPORT_RADIUS)
    return weight * excess**2

deg_lo, deg_hi = bbd["Photodeg_pct"].min(), bbd["Photodeg_pct"].max()
eeo_lo, eeo_hi = bbd["EEO_kWh_m3_order"].min(), bbd["EEO_kWh_m3_order"].max()
sec_lo, sec_hi = bbd["SEC_total_Wh_mg"].min(), bbd["SEC_total_Wh_mg"].max()

def desirabilities(g, c, t):
    deg = rsm_predict("Photodeg_pct", g, c, t)
    eeo = rsm_predict("EEO_kWh_m3_order", g, c, t)
    sec = rsm_predict("SEC_total_Wh_mg", g, c, t)
    d_deg = np.clip((deg - deg_lo) / (deg_hi - deg_lo), 1e-6, 1)
    d_eeo = np.clip((eeo_hi - eeo) / (eeo_hi - eeo_lo), 1e-6, 1)
    d_sec = np.clip((sec_hi - sec) / (sec_hi - sec_lo), 1e-6, 1)
    return d_deg, d_eeo, d_sec

def desirability_obj(params):
    g, c, t = params
    d_deg, d_eeo, d_sec = desirabilities(g, c, t)
    D = (d_deg * d_eeo * d_sec) ** (1/3)
    return -D + support_penalty(g, c, t)

res_desire = differential_evolution(desirability_obj, BOUNDS_REAL, seed=RNG, popsize=15, maxiter=300, tol=1e-10)
g_d, c_d, t_d = res_desire.x

def weighted_sum_obj(params, w=(1/3,1/3,1/3)):
    g, c, t = params
    d_deg, d_eeo, d_sec = desirabilities(g, c, t)
    return -(w[0]*d_deg + w[1]*d_eeo + w[2]*d_sec) + support_penalty(g, c, t)

res_ws = differential_evolution(weighted_sum_obj, BOUNDS_REAL, seed=RNG, popsize=15, maxiter=300, tol=1e-10)
g_w, c_w, t_w = res_ws.x

def tchebycheff_obj(params, ideal):
    g, c, t = params
    d_deg, d_eeo, d_sec = desirabilities(g, c, t)
    id_deg, id_eeo, id_sec = ideal
    return max(abs(id_deg-d_deg), abs(id_eeo-d_eeo), abs(id_sec-d_sec)) + support_penalty(g, c, t)

res_tcheb = differential_evolution(tchebycheff_obj, BOUNDS_REAL, args=((1,1,1),), seed=RNG, popsize=15, maxiter=300, tol=1e-10)
g_t, c_t, t_t = res_tcheb.x

for name, (g,c,t) in [("Desirability (DE)", (g_d,c_d,t_d)), ("Weighted sum (equal)", (g_w,c_w,t_w)),
                       ("Tchebycheff", (g_t,c_t,t_t))]:
    print(f"{name}: Glucose={g:.3f} CuSO4={c:.4f} Time={t:.1f}  ||X||={coded_norm(g,c,t):.3f} -> "
          f"Deg={rsm_predict('Photodeg_pct',g,c,t):.2f}%  EEO={rsm_predict('EEO_kWh_m3_order',g,c,t):.1f}  "
          f"SEC={rsm_predict('SEC_total_Wh_mg',g,c,t):.1f}")

try:
    from pymoo.core.problem import Problem
    from pymoo.algorithms.moo.nsga2 import NSGA2
    from pymoo.optimize import minimize as pymoo_minimize

    class TiO2Problem3(Problem):
        def __init__(self):
            super().__init__(n_var=3, n_obj=3, n_constr=1,
                              xl=np.array([0.5,0.02,120]), xu=np.array([2.0,0.08,180]))
        def _evaluate(self, X, out, *a, **k):
            g,c,t = X[:,0], X[:,1], X[:,2]
            f1 = -rsm_predict("Photodeg_pct", g,c,t)
            f2 = rsm_predict("EEO_kWh_m3_order", g,c,t)
            f3 = rsm_predict("SEC_total_Wh_mg", g,c,t)
            g1 = coded_norm(g,c,t) - BBD_SUPPORT_RADIUS
            out["F"] = np.column_stack([f1,f2,f3])
            out["G"] = np.column_stack([g1])

    algorithm = NSGA2(pop_size=100)
    result = pymoo_minimize(TiO2Problem3(), algorithm, ("n_gen", 150), seed=RNG, verbose=False)
    pareto_X, pareto_F = result.X, result.F
    pareto_df = pd.DataFrame(pareto_X, columns=["Glucose_g","CuSO4_g","Time_min"])
    pareto_df["Photodeg_pct"] = -pareto_F[:,0]
    pareto_df["EEO_kWh_m3_order"] = pareto_F[:,1]
    pareto_df["SEC_total_Wh_mg"] = pareto_F[:,2]
    best_deg_row = pareto_df.loc[pareto_df["Photodeg_pct"].idxmax()]
    best_eeo_row = pareto_df.loc[pareto_df["EEO_kWh_m3_order"].idxmin()]
    print(f"\nNSGA-II Pareto front: {len(pareto_df)} solutions")
    nsga2_available = True
except Exception as e:
    print(f"[NSGA-II skipped: {e}]")
    nsga2_available = False

opt_summary = pd.DataFrame([
    {"Method": "Desirability (DE)", "Glucose_g": g_d, "CuSO4_g": c_d, "Time_min": t_d,
     "Photodeg_pct": rsm_predict("Photodeg_pct",g_d,c_d,t_d), "EEO_kWh_m3_order": rsm_predict("EEO_kWh_m3_order",g_d,c_d,t_d),
     "SEC_total_Wh_mg": rsm_predict("SEC_total_Wh_mg",g_d,c_d,t_d)},
    {"Method": "Weighted sum (equal)", "Glucose_g": g_w, "CuSO4_g": c_w, "Time_min": t_w,
     "Photodeg_pct": rsm_predict("Photodeg_pct",g_w,c_w,t_w), "EEO_kWh_m3_order": rsm_predict("EEO_kWh_m3_order",g_w,c_w,t_w),
     "SEC_total_Wh_mg": rsm_predict("SEC_total_Wh_mg",g_w,c_w,t_w)},
    {"Method": "Tchebycheff", "Glucose_g": g_t, "CuSO4_g": c_t, "Time_min": t_t,
     "Photodeg_pct": rsm_predict("Photodeg_pct",g_t,c_t,t_t), "EEO_kWh_m3_order": rsm_predict("EEO_kWh_m3_order",g_t,c_t,t_t),
     "SEC_total_Wh_mg": rsm_predict("SEC_total_Wh_mg",g_t,c_t,t_t)},
])
if nsga2_available:
    opt_summary = pd.concat([opt_summary, pd.DataFrame([
        {"Method": "NSGA-II (best degradation)", **best_deg_row.to_dict()},
        {"Method": "NSGA-II (best EEO)", **best_eeo_row.to_dict()},
    ])], ignore_index=True)
print("\n", opt_summary.round(3).to_string(index=False))

# Also predict at the actually-tested new condition, for the validation section
NEW_COND = (1.57, 0.072, 140)
pred_at_new = {resp: rsm_predict(resp, *NEW_COND) for resp in RESPONSES}
print("\nModel prediction at the tested new condition (1.57/0.072/140):", pred_at_new)

# ---------------------------------------------------------------------------
# 6. MECHANISTIC ABLATION KINETICS (K1-K6)
# ---------------------------------------------------------------------------

print("\n" + "=" * 78 + "\nSECTION 6: MECHANISTIC ABLATION\n" + "=" * 78)
TIME_POINTS = np.array([0,15,30,45,60], dtype=float)

def pseudo_first_order_fit(a0, a_series):
    y = np.log(a0 / np.array(a_series))
    X = sm.add_constant(TIME_POINTS)
    model = sm.OLS(y, X).fit()
    return model.params[1], model.rsquared

kin_rows = []
for cond, g in KINETICS_RAW.groupby("Condition", sort=False):
    method = g["Method"].iloc[0]
    ks, r2s, removals = [], [], []
    for _, row in g.iterrows():
        series = [row["A0"], row["A15"], row["A30"], row["A45"], row["A60"]]
        k_reg, r2 = pseudo_first_order_fit(row["A0"], series)
        removals.append(photodegradation_pct(row["A0"], row["A60"]))
        ks.append(k_reg); r2s.append(r2)
    kin_rows.append({
        "Condition": cond, "Method": method,
        "Photodeg_60min_pct_mean": np.mean(removals), "Photodeg_60min_pct_sd": np.std(removals, ddof=1),
        "k_regression_min-1_mean": np.mean(ks), "k_regression_min-1_sd": np.std(ks, ddof=1),
        "R2_mean": np.mean(r2s),
        "EEO_mean": np.mean([eeo_kwh_m3_order(LED_POWER_W, PHOTO_TIME_MIN, SAMPLE_VOLUME_L, row["A0"], row["A60"])
                              for _, row in g.iterrows()]),
    })
kin_summary = pd.DataFrame(kin_rows)
print(kin_summary.round(5).to_string(index=False))

# ---------------------------------------------------------------------------
# 7. EXPERIMENTAL VALIDATION (V1-V3) vs NEW-MODEL PREDICTION & DESIGN CENTER
# ---------------------------------------------------------------------------

print("\n" + "=" * 78 + "\nSECTION 7: VALIDATION (V1-V3)\n" + "=" * 78)
center_bbd = bbd[bbd["Run"].isin([13,14,15])]
center_photodeg = center_bbd["Photodeg_pct"].values
center_kapp = center_bbd["k_app_min-1"].values

val = build_response_table(VALIDATION_RAW)

# Pooled direct replicates of the new optimum (V1-V3), computed first so V4-V7 can be
# expressed as a recovery relative to it (the actual sensitivity question for a
# single-factor perturbation set: how much does moving off the optimum cost?).
pooled = val[val["Set"].isin(["V1", "V2", "V3"])]
pooled_deg = pooled["Photodeg_pct"].values
t_pool_c, p_pool_c = stats.ttest_ind(center_photodeg, pooled_deg, equal_var=False)
pooled_recovery_vs_center = pooled_deg.mean() / center_photodeg.mean() * 100
pooled_recovery_vs_model = pooled_deg.mean() / pred_at_new["Photodeg_pct"] * 100

set_conditions = VALIDATION_RAW.groupby("Set", sort=False)[["Glucose_g", "CuSO4_g", "Time_min"]].first()

val_rows = []
for vset, g in val.groupby("Set", sort=False):
    photodeg, kapp = g["Photodeg_pct"].values, g["k_app_min-1"].values
    t_deg_c, p_deg_c = stats.ttest_ind(center_photodeg, photodeg, equal_var=False)
    t_k_c, p_k_c = stats.ttest_ind(center_kapp, kapp, equal_var=False)
    cnd = set_conditions.loc[vset]
    pred_at_this_point = rsm_predict("Photodeg_pct", cnd["Glucose_g"], cnd["CuSO4_g"], cnd["Time_min"])
    val_rows.append({
        "Set": vset, "Label": VALIDATION_SET_LABELS[vset],
        "Glucose_g": float(cnd["Glucose_g"]), "CuSO4_g": float(cnd["CuSO4_g"]), "Time_min": float(cnd["Time_min"]),
        "Photodeg_mean": photodeg.mean(), "Photodeg_sd": photodeg.std(ddof=1),
        "Recovery_vs_center_pct": photodeg.mean() / center_photodeg.mean() * 100,
        "Recovery_vs_new_optimum_pct": photodeg.mean() / pooled_deg.mean() * 100,
        "Model_prediction_at_point_pct": pred_at_this_point,
        "Recovery_vs_model_pct": photodeg.mean() / pred_at_this_point * 100,
        "Welch_t_vs_center": t_deg_c, "p_vs_center": p_deg_c,
        "kapp_mean": kapp.mean(), "kapp_sd": kapp.std(ddof=1),
        "EEO_mean": g["EEO_kWh_m3_order"].mean(), "EEO_sd": g["EEO_kWh_m3_order"].std(ddof=1),
        "SEC_total_mean": g["SEC_total_Wh_mg"].mean(), "SEC_total_sd": g["SEC_total_Wh_mg"].std(ddof=1),
    })
val_summary = pd.DataFrame(val_rows)
print(f"\nDesign-center (Runs 13-15, n=3): Photodeg={center_photodeg.mean():.3f}+/-{center_photodeg.std(ddof=1):.3f}%")
print(f"Model prediction at new optimum (1.57/0.072/140): Photodeg={pred_at_new['Photodeg_pct']:.3f}%")
print(val_summary.round(4).to_string(index=False))

print(f"\nPooled V1-V3 (n=9): Photodeg={pooled_deg.mean():.3f}+/-{pooled_deg.std(ddof=1):.3f}%, "
      f"recovery vs center={pooled_recovery_vs_center:.1f}%, recovery vs model={pooled_recovery_vs_model:.1f}%, "
      f"Welch t (vs center)={t_pool_c:.3f}, p={p_pool_c:.4f}")

# ---------------------------------------------------------------------------
# 8. EXPORT ../data/manuscript_data.json (only when run as a script, not on import by tio2_figures.py)
# ---------------------------------------------------------------------------

DATA = "../data/manuscript_data.json"
LEGACY_KEYS = ("lca_midpoint", "lca_sensitivity", "lca_montecarlo", "table9_endpoints", "lca_recipe_correction",
               "lca_normalized", "energy_intensity", "lca_single_process")   # superseded LCA blocks; the LCA is tio2_lca_single_process.py
if __name__ == "__main__":
  OLD = json.load(open(DATA, encoding="utf-8"))
  merged = {k: v for k, v in OLD.items() if k not in LEGACY_KEYS}  # keep characterization / antibacterial / optimization-hull sections
  merged.update({
      "bbd": json.loads(bbd.drop(columns=["X1","X2","X3"]).to_json(orient="records")),
      "rsm_summary": json.loads(rsm_summary.to_json(orient="records")),
      "ml_summary": json.loads(ml_summary.to_json(orient="records")),
      "opt_summary": json.loads(opt_summary.to_json(orient="records")),
      "kin_summary": json.loads(kin_summary.to_json(orient="records")),
      "val_summary": json.loads(val_summary.to_json(orient="records")),
      "rsm_coefs": rsm_coefs,
      "anova": anova_dict,
      "center_photodeg_mean": float(center_photodeg.mean()), "center_photodeg_sd": float(center_photodeg.std(ddof=1)),
      "center_kapp_mean": float(center_kapp.mean()), "center_kapp_sd": float(center_kapp.std(ddof=1)),
      "pareto": pareto_df.to_dict(orient="records") if nsga2_available else [],
      "bbd_raw": json.loads(BBD_RAW.to_json(orient="records")),
      "validation_raw": json.loads(VALIDATION_RAW.to_json(orient="records")),
      "kinetics_raw": json.loads(KINETICS_RAW.to_json(orient="records")),
      "new_optimum_condition": {"Glucose_g": NEW_COND[0], "CuSO4_g": NEW_COND[1], "Time_min": NEW_COND[2]},
      "new_optimum_model_prediction": pred_at_new,
      "pooled_v1v3_recovery_vs_center_pct": pooled_recovery_vs_center,
      "pooled_v1v3_recovery_vs_model_pct": pooled_recovery_vs_model,
      "pooled_v1v3_welch_t_vs_center": float(t_pool_c), "pooled_v1v3_p_vs_center": float(p_pool_c),
      "ppfd_umol_m2_s": PPFD_UMOL_M2_S,
  })

  with open(DATA, "w", encoding="utf-8") as f:
      json.dump(merged, f, indent=1, default=float)
  print("\nSaved ../data/manuscript_data.json (run tio2_reoptimize_hull.py next for the convex-hull optimization)")
