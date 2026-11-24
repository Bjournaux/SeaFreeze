"""
Global fit of the revised Redlich-Kister H2O-NH3 model to ice liquidus data.

Run from anywhere:   python run_fit.py
Reads the literature CSVs in ../For_Baptiste/data/literature (not modified).
Writes everything to ./outputs.
"""

import json
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nh3_rk as m

LIT = HERE.parent / "For_Baptiste" / "data" / "literature"
ORIG = HERE.parent / "For_Baptiste"
OUT = HERE / "outputs"
OUT.mkdir(exist_ok=True)

# ----------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------

LELIWA_UNITS = "wt"      # NH3 wt%, confirmed by the paper ("mol" kept only for the sensitivity test)
BOUND = 500.0            # |parameter| bound; the run aborts if one is reached
# Which ice each data point is fitted with:
#   "stable": high-pressure points (labelled III or V) use the ice with the lowest Gibbs
#             energy at the measured (P, T) among II, III, V, VI: ice II replaces ice III
#             at low temperature, ice VI replaces ice V at high pressure.
#             Ice Ih labels are kept: the 2 wt% points above 210 MPa follow the metastable
#             Ih melting curve (T still falling with P), not the ice III liquidus.
#   "label" : the ice named in the data file
ICE_ASSIGNMENT = "stable"

# Deliberate relabelling of published phase labels (a decision, not a rule): the series below
# are labelled ice III in the source but are treated as ice II.
#   Decision : B. Journaux, 2026-11-20 - "I think it was an error from the experimentalists."
#   Evidence : as ice III the 20 wt% series misses by -4 to -7 K with a pressure trend, as ice II by
#              +1.6, -0.3, -0.4 K. At these (P, T) ice II has the lower Gibbs energy. Choukroun &
#              Grasset (2010) themselves conclude that ice III cannot exist above 5-10 wt% NH3.
#   To document in the paper. Every label change is written to outputs/ice_label_changes.csv.
RELABEL = {
    "Chroukroun_x02_III": ("II", "relabelled III -> II by decision (B. Journaux, 2026-11-20): "
                                 "suspected misidentification in the source; 20 wt%"),
}

# Series whose phase is decided point by point: each point takes, among the candidate ices, the
# one giving the smallest |T(model) - T(observed)|.  The labels depend on the fit and the fit on the
# labels, so main() repeats until they stop changing.
#   Decision : B. Journaux, 2026-11-20, for the Choukroun & Grasset 10 wt% "ice III" series, which
#              lies near the II-III boundary and fits either ice within ~1 K.
#   Caveat   : this selects labels on the outcome, so the misfit of these points is a lower bound,
#              and ice II is the lower-Gibbs-energy ice at all four points (ice III there would
#              be metastable).
BEST_LABEL = {
    "Chroukroun_x01_III": ("II", "III"),
}
_BEST_STATE = {name: None for name in BEST_LABEL}     # per-point ices chosen so far (None: first candidate)
BEST_NOTE = ("best-fitting ice among {cands} for this point, by decision (B. Journaux, 2026-11-20); "
             "source label III; 10 wt%")
N_BOOT = 500
RNG = np.random.default_rng(0)

# Model ladder: parameter subsets fitted globally to all ices
LADDER = {
    "A1  a0":                 ("a0",),
    "A2  a0 a1":              ("a0", "a1"),
    "B1  a0 + v0":            ("a0", "v0"),
    "C1  a0 + b0":            ("a0", "b0"),
    "D1  a0 + b0 v0":         ("a0", "b0", "v0"),
    "A3  a0 a1 a2":           ("a0", "a1", "a2"),
    "B   a0 a1 + v0":         ("a0", "a1", "v0"),
    "C   a0 a1 + b0":         ("a0", "a1", "b0"),
    "D   a0 a1 + b0 v0":      ("a0", "a1", "b0", "v0"),
    "D3  a0 a1 a2 + b0 v0":   ("a0", "a1", "a2", "b0", "v0"),
    "E   a0 a1 + b0 v0 v1":   ("a0", "a1", "b0", "v0", "v1"),
    "F   a0 a1 + b0 b1 v0":   ("a0", "a1", "b0", "b1", "v0"),
    "G   a0 a1 + b0 b1 v0 v1": ("a0", "a1", "b0", "b1", "v0", "v1"),
}


LELIWA_FIT = "Leliwa-Kopystynski 2002 (Table III fit)"


def make_datasets(units=LELIWA_UNITS):
    def leliwa_x(number):
        return float(m.wt_to_x(number)) if units == "wt" else number / 100.0

    # (file, ice, x_NH3, source, used in fit, note)
    return [
        ("amm2pct_ih.csv",        "Ih",  float(m.wt_to_x(2)),   "Choukroun & Grasset 2010", True,  "2 wt%"),
        ("amm8pct_ih.csv",        "Ih",  float(m.wt_to_x(8)),   "Choukroun & Grasset 2010", True,  "8 wt%"),
        ("amm11pct_ih.csv",       "Ih",  float(m.wt_to_x(11)),  "Choukroun & Grasset 2010", True,  "11 wt%"),
        ("Chroukroun_x01_Ih.csv", "Ih",  float(m.wt_to_x(10)),  "Choukroun & Grasset 2010", True,  "10 wt%"),
        ("Chroukroun_x02_Ih.csv", "Ih",  float(m.wt_to_x(20)),  "Choukroun & Grasset 2010", True,  "20 wt%"),
        ("Leliwa_x3.34.csv",      "Ih",  leliwa_x(3.34),        "Leliwa-Kopystynski 2002",  True,  "3.34"),
        ("Leliwa_x4.72.csv",      "Ih",  leliwa_x(4.72),        "Leliwa-Kopystynski 2002",  True,  "4.72"),
        ("Leliwa_x17.6.csv",      "Ih",  leliwa_x(17.6),        "Leliwa-Kopystynski 2002",  True,  "17.6"),
        ("Leliwa_x22.9.csv",      "Ih",  leliwa_x(22.9),        "Leliwa-Kopystynski 2002",  True,  "22.9"),
        # Series missing from the digitised folder, generated from the paper's own fits (Table III).
        # Same number of points as a digitised series, so they do not outweigh the raw data.
        ("data_added/Leliwa_x11.1_Ih_TableIII.csv",  "Ih",  leliwa_x(11.1), LELIWA_FIT, True, "11.1, Table III fit"),
        ("data_added/Leliwa_x3.34_III_TableIII.csv", "III", leliwa_x(3.34), LELIWA_FIT, True, "3.34 ice III, Table III fit"),
        ("data_added/Leliwa_x4.72_III_TableIII.csv", "III", leliwa_x(4.72), LELIWA_FIT, True, "4.72 ice III, Table III fit"),
        # amm014: 10-12 K below the pure-water ice III melting curve, impossible for x = 0.026.
        # Composition (or phase) label is wrong -> shown, not fitted.
        ("amm014_III.csv",        "III", 0.0264,                "Choukroun & Grasset 2010", False, "x as given in original main.py; inconsistent"),
        ("amm015_III.csv",        "III", 0.0741,                "Choukroun & Grasset 2010", True,  "x as given in original main.py"),
        ("Chroukroun_x01_III.csv", "III", float(m.wt_to_x(10)), "Choukroun & Grasset 2010", True,  "10 wt%"),
        ("Chroukroun_x02_III.csv", "III", float(m.wt_to_x(20)), "Choukroun & Grasset 2010", True,  "20 wt%"),
        ("amm015_V.csv",          "V",   0.0741,                "Choukroun & Grasset 2010", True,  "x as given in original main.py"),
        ("Chroukroun_x01_V.csv",  "V",   float(m.wt_to_x(10)),  "Choukroun & Grasset 2010", True,  "10 wt%"),
        ("Chroukroun_x02_V.csv",  "V",   float(m.wt_to_x(20)),  "Choukroun & Grasset 2010", True,  "20 wt%"),
        # Not loaded:
        #  Leliwa_x11.1.csv  is the 4.72 wt% data (byte-identical copy of Leliwa_x4.72.csv)
        #  amm013_II*.csv    no composition given, two inconsistent digitisations
    ]


# ----------------------------------------------------------------------------
# Data
# ----------------------------------------------------------------------------

def load_data(units=LELIWA_UNITS, verbose=True):
    rows = []
    for i, (fname, ice, x, source, used, note) in enumerate(make_datasets(units)):
        path = (HERE / fname) if fname.startswith("data_added/") else (LIT / fname)
        d = np.array([[float(v) for v in line.split(",")] for line in path.read_text().splitlines()
                      if line.strip() and line.lstrip()[0] in "0123456789.-"])
        for P, T in d:
            rows.append(dict(set=Path(fname).name.replace(".csv", ""), iset=i, ice=ice, x=x,
                             source=source, used=used, P=P, T=T))
    df = pd.DataFrame(rows)
    dG = np.full(len(df), np.nan); dS = np.full(len(df), np.nan)
    for ice in df.ice.unique():
        k = (df.ice == ice).values
        dG[k], dS[k] = m.fusion(df.P.values[k], df["T"].values[k], ice)
    df["dG"], df["dS"] = dG, dS
    df["ice_label"] = df.ice                      # label in the source, never modified
    df["label_change"] = ""
    # 1. deliberate relabelling (applies whatever ICE_ASSIGNMENT is)
    fixed = np.zeros(len(df), bool)
    for name, (new_ice, why) in RELABEL.items():
        k = (df.set == name).values
        g, e = m.fusion(df.P.values[k], df["T"].values[k], new_ice)
        df.loc[k, ["dG", "dS"]] = np.column_stack([g, e])
        df.loc[k, "ice"] = new_ice
        df.loc[k, "label_change"] = why
        fixed |= k
    # 1b. series decided point by point (labels come from the previous pass of main())
    for name, cands in BEST_LABEL.items():
        idx = np.where((df.set == name).values)[0]
        chosen = _BEST_STATE[name] or [cands[0]] * len(idx)
        for i, ice in zip(idx, chosen):
            g, e = m.fusion(df.P.values[[i]], df["T"].values[[i]], ice)
            df.loc[i, ["dG", "dS"]] = [g[0], e[0]]
            df.loc[i, "ice"] = ice
            df.loc[i, "label_change"] = BEST_NOTE.format(cands="/".join(cands)) + \
                ("" if ice != df.ice_label[i] else " -> source label kept")
        fixed[idx] = True
    # 2. stable-ice rule for the remaining high-pressure points
    if ICE_ASSIGNMENT == "stable":
        hp = (df.ice_label != "Ih").values & ~fixed
        for ice in ("II", "III", "V", "VI"):
            g, e = m.fusion(df.P.values, df["T"].values, ice)
            better = hp & np.isfinite(g) & ~(g >= df.dG.values)   # lower G_ice, or label was NaN
            df.loc[better, "dG"] = g[better]
            df.loc[better, "dS"] = e[better]
            df.loc[better, "ice"] = ice
    df["best_label_series"] = df.set.isin(list(BEST_LABEL))
    rule = (df.ice != df.ice_label).values & ~fixed
    df.loc[rule, "label_change"] = ("stable-ice rule: " + df.ice_label[rule] + " -> " + df.ice[rule]
                                    + " (lower Gibbs energy at the measured P, T)")
    moved = df[df.ice != df.ice_label]
    if verbose and len(moved):
        print("Points fitted with an ice other than the one labelled in the source:")
        print(moved.groupby(["set", "ice_label", "ice", "label_change"]).size().to_string())
    bad = ~np.isfinite(df.dG)
    if bad.any() and verbose:
        print("Points outside the SeaFreeze domain of their ice (dropped, not clipped):")
        print(df.loc[bad, ["set", "ice", "P", "T"]].to_string(index=False))
    return df[~bad].reset_index(drop=True)


# ----------------------------------------------------------------------------
# Fit (the residual is linear in the parameters)
# ----------------------------------------------------------------------------

def linear_system(df, names):
    """Residual in kelvin:  r = A p - y  ~  T_model - T_obs."""
    RT = m.R * df["T"].values
    A = m.design_matrix(df.P.values, df["T"].values, df.x.values, names)
    A = -(RT / df.dS.values)[:, None] * A
    y = -(df.dG.values - RT * np.log1p(-df.x.values)) / df.dS.values
    return A, y


def fit(df, names, robust=True):
    A, y = linear_system(df, names)
    p0 = np.linalg.lstsq(A, y, rcond=None)[0]
    p0 = np.clip(p0, -0.99 * BOUND, 0.99 * BOUND)
    res = least_squares(lambda p: A @ p - y, p0, jac=lambda p: A,
                        loss="soft_l1" if robust else "linear", f_scale=1.0,
                        bounds=(-BOUND, BOUND))
    return res.x, A @ res.x - y, A


def rms(r):
    return float(np.sqrt(np.mean(np.asarray(r)**2)))


def leave_one_set_out(df, names):
    r_all = []
    for s in df.set.unique():
        train, test = df[df.set != s], df[df.set == s]
        p, _, _ = fit(train, names)
        A, y = linear_system(test, names)
        r_all.append(A @ p - y)
    return rms(np.concatenate(r_all))


def bootstrap(df, names, n=N_BOOT):
    sets = df.set.unique()
    P = []
    for _ in range(n):
        pick = RNG.choice(sets, size=len(sets), replace=True)
        d = pd.concat([df[df.set == s] for s in pick], ignore_index=True)
        A, _ = linear_system(d, names)
        if np.linalg.matrix_rank(A) < len(names):
            continue
        P.append(fit(d, names)[0])
    return np.array(P)


# ----------------------------------------------------------------------------
# Original model (For_Baptiste) evaluated for comparison
# ----------------------------------------------------------------------------

def original_model_dT(df):
    """T_model - T_obs of the original code (its own tables and fit). NaN on failure."""
    try:
        cp = types.ModuleType("CoolProp"); cpp = types.ModuleType("CoolProp.CoolProp")
        cpp.PropsSI = lambda *a, **k: float("nan"); cp.CoolProp = cpp
        sys.modules.setdefault("CoolProp", cp); sys.modules.setdefault("CoolProp.CoolProp", cpp)
        import importlib.util, os, warnings
        cwd = os.getcwd(); os.chdir(HERE)           # utils creates ./outputs at import
        spec = importlib.util.spec_from_file_location("orig_utils", ORIG / "utils.py")
        u = importlib.util.module_from_spec(spec)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            spec.loader.exec_module(u)
            sfd = ORIG / "data" / "seafreeze"
            suffix = {"Ih": "", "II": "_II", "III": "_III", "V": "_V", "VI": "_VI"}
            tabs = {i: u.initialize_fusion_model(str(sfd / f"GibbsEnergywater_IAPWS95_kg{s}.txt"),
                                                 str(sfd / f"GibbsEnergy{i}_kg.txt"))
                    for i, s in suffix.items()}
            orig_sets = [("amm2pct_ih.csv", 0.0211), ("amm8pct_ih.csv", 0.0842),
                         ("amm11pct_ih.csv", 0.1156), ("Leliwa_x3.34.csv", 0.0334),
                         ("Leliwa_x4.72.csv", 0.0472), ("Leliwa_x11.1.csv", 0.111),
                         ("Leliwa_x17.6.csv", 0.1843)]
            Pa, Ta, xa = u.load_liquidus_datasets([(str(LIT / f), x) for f, x in orig_sets])
            L = u.fit_margules(Pa, Ta, xa, *tabs["Ih"]).x
            dT = [u.find_liquidus_temperature(r.P, r.x, *tabs[r.ice], L_fit=L) - r["T"]
                  for _, r in df.iterrows()]
        os.chdir(cwd)
        for k in ("CoolProp", "CoolProp.CoolProp"):
            if sys.modules.get(k) in (cp, cpp):
                del sys.modules[k]
        return np.array(dT), L
    except Exception as e:                           # comparison only
        print("Original-model comparison skipped:", repr(e))
        return np.full(len(df), np.nan), None


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main(_pass=1):
    df = load_data(verbose=(_pass == 1))
    dfit = df[df.used].reset_index(drop=True)
    n = len(dfit)
    print(f"\n{n} points in {dfit.set.nunique()} series; Leliwa compositions read as {LELIWA_UNITS}%")
    print(f"x range {dfit.x.min():.3f}-{dfit.x.max():.3f}, P {dfit.P.min():.0f}-{dfit.P.max():.0f} MPa, "
          f"T {dfit['T'].min():.1f}-{dfit['T'].max():.1f} K")

    # ---- ideal solution baseline
    _, y = linear_system(dfit, ("a0",))
    print(f"Ideal solution (no RK): RMS = {rms(-y):.2f} K")

    # ---- model ladder
    ladder = []
    for label, names in LADDER.items():
        p, r, A = fit(dfit, names)
        k = len(names)
        sse = float(np.sum(r**2))
        ladder.append(dict(model=label, k=k, rms_K=rms(r),
                           AIC=n * np.log(sse / n) + 2 * k,
                           BIC=n * np.log(sse / n) + k * np.log(n),
                           loso_rms_K=leave_one_set_out(dfit, names),
                           on_bound=bool(np.any(np.abs(p) > 0.98 * BOUND)),
                           cond=float(np.linalg.cond(A / np.linalg.norm(A, axis=0)))))
    ladder = pd.DataFrame(ladder)
    print("\nModel ladder (RMS in K; loso = leave-one-series-out prediction):")
    print(ladder.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    ladder.to_csv(OUT / "model_ladder.csv", index=False)

    # ---- selection: smallest model within 5 % of the best leave-one-series-out RMS
    ok = ladder[~ladder.on_bound]
    best = ok.loso_rms_K.min()
    chosen = ok[ok.loso_rms_K <= 1.05 * best].sort_values(["k", "loso_rms_K"]).iloc[0].model
    names = LADDER[chosen]
    print(f"\nSelected model: {chosen}")

    p_red, r, A = fit(dfit, names)
    if np.any(np.abs(p_red) > 0.98 * BOUND):
        raise RuntimeError("A parameter is on its bound: the fit is not constrained.")
    p = m.full_params(p_red, names)

    boot = bootstrap(dfit, names)
    sig = boot.std(axis=0)
    lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
    corr = np.corrcoef(boot.T) if len(names) > 1 else np.array([[1.0]])
    print("\nParameters (bootstrap over series):")
    for nme, v, s, a, b in zip(names, p_red, sig, lo, hi):
        print(f"  {nme} = {v:+9.3f} +/- {s:6.3f}   95% [{a:+.3f}, {b:+.3f}]")
    print("Correlation matrix:\n", pd.DataFrame(corr, index=names, columns=names).round(2).to_string())

    # ---- liquidus tables, true temperature misfit per series
    tables = m.FusionTables()

    # ---- point-by-point labels: keep the best-fitting candidate, repeat until stable
    changed = False
    for name, cands in BEST_LABEL.items():
        rows = df[df.set == name]
        best = []
        for _, row in rows.iterrows():
            dT = {ice: abs(m.liquidus_T(row.P, row.x, ice, p, tables) - row["T"]) for ice in cands}
            best.append(min(dT, key=lambda k: dT[k] if np.isfinite(dT[k]) else np.inf))
        if best != list(rows.ice):
            _BEST_STATE[name] = best
            changed = True
        print(f"\nBest-fitting labels for {name} (pass {_pass}): {best}")
    if changed and _pass < 8:
        print("Labels changed: refitting.\n")
        return main(_pass + 1)
    if changed:
        raise RuntimeError("point-by-point labels did not converge")
    df["dT_revised"] = [m.liquidus_T(r.P, r.x, r.ice, p, tables) - r["T"] for _, r in df.iterrows()]
    df["dT_revised_label"] = [m.liquidus_T(r.P, r.x, r.ice_label, p, tables) - r["T"] for _, r in df.iterrows()]
    df["dT_ideal"] = [m.liquidus_T(r.P, r.x, r.ice, np.zeros(9), tables) - r["T"] for _, r in df.iterrows()]
    df["dT_original"], L_orig = original_model_dT(df)
    stab = [m.liquidus_stable(r.P, r.x, p, tables) for _, r in df.iterrows()]
    df["ice_stable_model"] = [s_[1] for s_ in stab]
    df["dT_revised_stable"] = [s_[0] for s_ in stab] - df["T"].values
    diff = df[df.ice != df.ice_stable_model]
    if len(diff):
        print("\nPoints whose labelled ice is not the model's stable liquidus phase "
              "(fitted with the labelled ice):")
        print(diff.groupby(["set", "ice", "ice_stable_model"]).size().to_string())

    # ---- sensitivity to the Leliwa composition units
    alt = "mol" if LELIWA_UNITS == "wt" else "wt"
    d_alt = load_data(alt, verbose=False); d_alt = d_alt[d_alt.used].reset_index(drop=True)
    _, r_alt, _ = fit(d_alt, names)
    lel = dfit.set.str.startswith("Leliwa").values
    print(f"\nLeliwa units test, same model: RMS all / Leliwa-only [K]: "
          f"{LELIWA_UNITS}% {rms(r):.2f} / {rms(r[lel]):.2f}   |   {alt}% {rms(r_alt):.2f} / {rms(r_alt[lel]):.2f}")

    g = df.groupby("set", sort=False)
    mis = pd.DataFrame({
        "ice": g.ice.agg(lambda v: "/".join(pd.unique(v))), "fitted": g.used.first(), "x_NH3": g.x.first().round(4),
        "wt_pct": g.x.first().apply(lambda v: round(float(m.x_to_wt(v)), 2)),
        "n": g.size(), "P_min": g.P.min().round(0), "P_max": g.P.max().round(0),
        "T_min": g["T"].min().round(1),
        "mean_dT_original": g.dT_original.mean().round(2),
        "mean_dT_ideal": g.dT_ideal.mean().round(2),
        "mean_dT_revised": g.dT_revised.mean().round(2),
        "max_abs_dT_revised": g.dT_revised.apply(lambda v: np.nanmax(np.abs(v))).round(2),
    })
    print("\nMisfit per series, T_model - T_obs [K]:")
    print(mis.to_string())
    mis.to_csv(OUT / "misfit_by_dataset.csv")
    du = df[df.used]
    print(f"\nRMS over fitted points [K]: original {rms(du.dT_original.dropna()):.2f} | "
          f"ideal {rms(du.dT_ideal.dropna()):.2f} | revised {rms(du.dT_revised.dropna()):.2f}")
    for ice in du.ice.unique():
        d = du[du.ice == ice]
        print(f"  ice {ice:3s}: original {rms(d.dT_original.dropna()):.2f} | revised {rms(d.dT_revised.dropna()):.2f}")

    # ---- implied excess properties (sanity only)
    print("\nImplied excess properties at 240 K, 0.1 MPa (not fitted to such data):")
    for x in (0.1, 0.2):
        Gex, Hex, Vex = m.excess_properties(m.P_REF, m.T_REF, x, p)
        print(f"  x={x}: G_ex={Gex:+.0f} J/mol  H_ex={Hex:+.0f} J/mol  V_ex={Vex:+.2f} cm3/mol  "
              f"ln gamma_w={float(m.ln_gamma_w(m.P_REF, m.T_REF, x, p)):+.4f}")

    # ---- save parameters
    import seafreeze
    meta = dict(
        model="Redlich-Kister, L_k = a_k + b_k (T_ref/T - 1) + v_k (P - P_ref)/1000",
        selected=chosen, T_ref_K=m.T_REF, P_ref_MPa=m.P_REF,
        parameters={nme: float(v) for nme, v in zip(m.PARAM_NAMES, p)},
        fitted=list(names),
        bootstrap_sigma={nme: float(s) for nme, s in zip(names, sig)},
        bootstrap_95pct={nme: [float(a), float(b)] for nme, a, b in zip(names, lo, hi)},
        rms_K=rms(df[df.used].dT_revised.dropna()),
        n_points=int(n), n_series=int(dfit.set.nunique()),
        x_range=[float(dfit.x.min()), float(dfit.x.max())],
        P_range_MPa=[float(dfit.P.min()), float(dfit.P.max())],
        T_range_K=[float(dfit["T"].min()), float(dfit["T"].max())],
        leliwa_units=LELIWA_UNITS,
        water_reference="SeaFreeze 'water1' for T >= 240 K; constant-Cp extrapolation of "
                        "'water1' (G, S, Cp at 240 K) below 240 K",
        ices="SeaFreeze Ih, II, III, V, VI",
        seafreeze_version=getattr(seafreeze, "__version__", "unknown"),
        original_model_L=None if L_orig is None else [float(v) for v in L_orig],
    )
    (OUT / "fit_parameters.json").write_text(json.dumps(meta, indent=2))
    df.to_csv(OUT / "data_with_residuals.csv", index=False)
    log = df[(df.ice != df.ice_label) | df.best_label_series][["set", "source", "x", "P", "T", "ice_label", "ice", "used",
                                      "label_change", "dT_revised_label", "dT_revised"]]
    log = log.rename(columns={"x": "x_NH3", "P": "P_MPa", "T": "T_K", "ice_label": "ice_in_source",
                              "ice": "ice_used", "used": "used_in_fit",
                              "dT_revised_label": "dT_K_with_source_ice", "dT_revised": "dT_K_with_ice_used"})
    log.round(4).to_csv(OUT / "ice_label_changes.csv", index=False)
    print(f"\n{int((df.ice != df.ice_label).sum())} points use an ice other than the source label "
          f"({len(log)} rows logged); see outputs/ice_label_changes.csv")

    make_figures(df, p, tables)
    export_liquidus(p, tables)
    return df, p


# ----------------------------------------------------------------------------
# Figures and exports
# ----------------------------------------------------------------------------

MARK = {"Choukroun & Grasset 2010": "o", "Leliwa-Kopystynski 2002": "P", LELIWA_FIT: "X"}


def make_figures(df, p, tables):
    Pg = np.linspace(0.1, 900, 361)
    cmap = plt.get_cmap("viridis_r")
    xmax = 0.25
    col = lambda x: cmap(0.08 + 0.9 * x / xmax)

    # --- phase diagram
    fig, ax = plt.subplots(figsize=(8.5, 6))
    ax.plot(Pg, [m.liquidus_stable(P, 0.0, p, tables)[0] for P in Pg], c="k", lw=2, label="pure water")
    for w in (2, 5, 8, 11, 15, 20, 23):
        x = float(m.wt_to_x(w))
        ax.plot(Pg, [m.liquidus_stable(P, x, p, tables)[0] for P in Pg], c=col(x), lw=1.6,
                label=f"{w} wt% (x = {x:.3f})")
    # stable pure-ice fields (lowest Gibbs energy), drawn below the pure-water melting curve
    Pf = np.linspace(0.1, 1000, 401); Tf = np.linspace(200, 300, 201)
    Gi = np.stack([np.where(np.isfinite(g), g, np.inf)
                   for g in (m._sf_grid(Pf, Tf, ice, ("G",))[0] for ice in m.ICES)])
    idx = np.argmin(Gi, axis=0).astype(float)
    Tm = np.array([m.liquidus_stable(P, 0.0, p, tables)[0] for P in Pf])
    idx[Tf[None, :] > Tm[:, None]] = np.nan
    for i in range(len(m.ICES)):
        ax.contour(Pf, Tf, (idx == i).astype(float).T, levels=[0.5], colors="gray",
                   linewidths=0.7, linestyles=":")
    for name, (px, ty) in {"Ih": (60, 205), "II": (330, 205), "III": (275, 249.5),
                           "V": (500, 258), "VI": (830, 205)}.items():
        ax.text(px, ty, f"ice {name}", color="gray", fontsize=8, ha="center")
    for (s, src), d in df.groupby(["set", "source"], sort=False):
        ax.scatter(d.P, d["T"], color=col(d.x.iloc[0]), marker=MARK[src], s=34,
                   edgecolor="k", linewidth=0.4, zorder=3)
    for src, mk in MARK.items():
        ax.scatter([], [], color="w", edgecolor="k", marker=mk, s=34, label=src)
    ax.set(xlabel="Pressure [MPa]", ylabel="Temperature [K]", xlim=(0, 1000), ylim=(200, 300),
           title="H$_2$O-NH$_3$ ice liquidus, revised Redlich-Kister model\n"
                 "(lines: stable liquidus; symbols coloured by composition; dotted: pure-ice fields)")
    ax.axhline(m.T_REF, c="gray", ls=":", lw=0.8)
    ax.text(990, m.T_REF + 0.5, "240 K: water reference extrapolated below", ha="right", fontsize=7, color="gray")
    ax.legend(frameon=False, fontsize=8, ncol=2, loc="upper left")
    ax.grid(True, ls=":", c="lightgray")
    fig.tight_layout(); fig.savefig(OUT / "phase_diagram_revised.pdf"); fig.savefig(OUT / "phase_diagram_revised.png", dpi=160)
    plt.close(fig)

    # --- residuals of the revised model, coloured by the ice each point is fitted with
    ICE_COL = {"Ih": "#0072B2", "II": "#D55E00", "III": "#009E73", "V": "#CC79A7", "VI": "#E69F00"}
    d = df[df.used]
    moved = d[d.ice != d.ice_label]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6), sharey=True)
    for ax, key, lab in zip(axes, ("P", "T", "x"), ("Pressure [MPa]", "Temperature [K]", "x(NH$_3$)")):
        ax.axhspan(-1, 1, color="0.93", zorder=0)
        ax.axhline(0, c="k", lw=0.8)
        # same point evaluated with the ice named in the data file (hollow), joined to the fitted one
        for _, r in moved.iterrows():
            ax.plot([r[key], r[key]], [r.dT_revised_label, r.dT_revised], c="0.6", lw=0.8, zorder=1)
            ax.scatter(r[key], r.dT_revised_label, marker=MARK[r.source], s=36, facecolor="white",
                       edgecolor=ICE_COL[r.ice_label], linewidth=1.2, zorder=2)
        for (ice, src), g in d.groupby(["ice", "source"]):
            ax.scatter(g[key], g.dT_revised, marker=MARK[src], s=36, color=ICE_COL[ice],
                       edgecolor="k", linewidth=0.3, zorder=3)
        ax.set_xlabel(lab); ax.grid(True, ls=":", c="lightgray")
    axes[0].set_ylabel("T$_{model}$ - T$_{obs}$ [K]")
    h = [plt.Line2D([], [], ls="", marker="o", color=c, markeredgecolor="k", markeredgewidth=0.3,
                    label=f"ice {i} (n = {(d.ice == i).sum()}, RMS {rms(d.dT_revised[d.ice == i]):.1f} K)")
         for i, c in ICE_COL.items() if (d.ice == i).any()]
    h += [plt.Line2D([], [], ls="", marker="o", markerfacecolor="white", markeredgecolor="0.3",
                     label="same point with the ice named in the file")]
    h += [plt.Line2D([], [], ls="", marker=mk, color="0.5", label=src) for src, mk in MARK.items()]
    fig.legend(handles=h, loc="upper center", ncol=4, frameon=False, fontsize=8)
    fig.suptitle("Residuals of the revised model by ice polymorph (grey band: +/- 1 K)", y=0.03, fontsize=9)
    fig.tight_layout(rect=(0, 0.04, 1, 0.88))
    fig.savefig(OUT / "residuals.pdf"); fig.savefig(OUT / "residuals.png", dpi=160)
    plt.close(fig)

    # --- water reference
    T = np.arange(200, 281, 1.0)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for P, c in zip((0.1, 200.0, 400.0, 600.0, 900.0), plt.get_cmap("plasma")(np.linspace(0, 0.85, 5))):
        G, _ = m.water_ref(P, T)
        Gi, Cpi = m._sf_scatter(P, T, "water_IAPWS95", ("G", "Cp"))
        axes[0].plot(T, Gi * m.M_H2O - G, c=c, label=f"{P:g} MPa")
        h = 0.5
        Ga, _ = m.water_ref(P, T + h); Gc, _ = m.water_ref(P, T - h)
        axes[1].plot(T, -T * (Ga - 2 * G + Gc) / h**2 / m.M_H2O, c=c)
        axes[1].plot(T, Cpi, c=c, ls="--")
    axes[0].set(xlabel="Temperature [K]", ylabel="G(IAPWS-95) - G(reference) [J/mol]",
                title="Original water tables vs. revised reference")
    axes[1].set(xlabel="Temperature [K]", ylabel="Cp [J/kg/K]", ylim=(-3000, 8000),
                title="Cp: revised reference (solid), IAPWS-95 (dashed)")
    for ax in axes:
        ax.axvline(m.T_REF, c="gray", ls=":"); ax.grid(True, ls=":", c="lightgray")
    axes[1].axhline(0, c="k", lw=0.8)
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "water_reference.pdf"); fig.savefig(OUT / "water_reference.png", dpi=160)
    plt.close(fig)


def export_liquidus(p, tables):
    Pg = np.concatenate([[0.1], np.arange(5.0, 900.1, 5.0)])
    for w in (0, 2, 5, 8, 10, 11, 15, 17.6, 20, 22.9):
        x = float(m.wt_to_x(w))
        rows = []
        for P in Pg:
            T, ice = m.liquidus_stable(P, x, p, tables)
            if not np.isfinite(T):
                rows.append(dict(P_MPa=P, T_liquidus_K=np.nan, ice=None)); continue
            lg = float(m.ln_gamma_w(P, T, x, p))
            rows.append(dict(P_MPa=P, T_liquidus_K=round(T, 3), ice=ice, x_NH3=round(x, 5), wt_pct=w,
                             ln_gamma_w=lg, a_w=float(np.exp(m.ln_a_w(P, T, x, p))),
                             mu_w_minus_G0w_J_mol=float(m.R * T * m.ln_a_w(P, T, x, p)),
                             water_ref_extrapolated=bool(T < m.T_REF)))
        pd.DataFrame(rows).to_csv(OUT / f"liquidus_{w:g}wt.csv", index=False)
    try:
        m.liquid_gibbs_total(100.0, 250.0, 0.1, p)
    except ImportError:
        print("\nCoolProp not installed: total G_liquid (needs pure NH3) not exported.")
    print(f"Outputs written to {OUT}")


if __name__ == "__main__":
    main()
