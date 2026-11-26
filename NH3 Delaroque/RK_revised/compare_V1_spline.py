"""
Compare the experimental 'NH3aq_RK' material (Redlich-Kister, fitted to ice liquidi) with the
NH3-H2O Gibbs spline already in the repository:

    Python/mlbspline/test/NH3_H2O_LBF.mat   variable spG_PTm_NH3_H2O_V1   ("V1")
    0-800 MPa, 239-401 K, 0-24.5 mol/kg; temperature axis is ln(T / 450 K)

    python "NH3 Delaroque/RK_revised/compare_V1_spline.py"      (from the repository root)
"""
import csv
import sys
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "Python"))
import seafreeze as sf                       # noqa: E402
import seafreeze.seafreeze as sfm            # noqa: E402
from mlbspline import load                   # noqa: E402
sys.path.insert(0, str(HERE))
from spline_eval import eval_bspline      # noqa: E402

warnings.simplefilter("ignore")
OUT = HERE / "outputs"
R, MW, M_NH3 = 8.314462618, sfm.mH2O_kgmol, sfm.mNH3_kgmol
NW = 1.0 / MW
wt2m = lambda w: (np.asarray(w, float) / (100.0 - np.asarray(w, float))) / M_NH3

V1 = load.loadSpline(str(REPO / "Python" / "mlbspline" / "test" / "NH3_H2O_LBF.mat"), "spG_PTm_NH3_H2O_V1")
TC = float(V1["Tc"])
M_V1 = float(V1["MW"])


def v1_props(P, T, m, mask=True):
    """rho, vel, Cp, muw from the V1 spline on a (nP, nT, nm) grid."""
    P = np.atleast_1d(np.asarray(P, float)); T = np.atleast_1d(np.asarray(T, float)); m = np.atleast_1d(np.asarray(m, float))
    x = np.empty(3, dtype=object); x[0], x[1], x[2] = P, np.log(T / TC), m
    d = lambda *o: eval_bspline(V1, x, list(o))     # NOT mlbspline's splev-based evaluator: see spline_eval.py
    Tm, mm = T[None, :, None], m[None, None, :]
    G, GP, GPP, Gt, Gtt, GPt, Gm = d(0, 0, 0), d(1, 0, 0), d(2, 0, 0), d(0, 1, 0), d(0, 2, 0), d(1, 1, 0), d(0, 0, 1)
    GTT, GPT = (Gtt - Gt) / Tm**2, GPt / Tm                     # chain rule for tau = ln(T/Tc)
    with np.errstate(all="ignore"):
        out = dict(rho=1e6 / GP, vel=np.sqrt(GP**2 / (GPT**2 / GTT - GPP)), Cp=-Tm * GTT,
                   muw=MW * (G - (1 + M_V1 * mm) * mm * Gm))
        # V1 oscillates outside the region it was fitted on (above ~500-600 MPa at low
        # temperature, and near its 239 K edge at low concentration): keep physical values only
        ok = (out["Cp"] > 2000) & (out["Cp"] < 6000) & (out["rho"] > 800) & (out["rho"] < 1400) \
            & (out["vel"] > 800) & (out["vel"] < 3500)
        out["physical"] = ok
        if mask:
            for k in ("rho", "vel", "Cp"):
                out[k] = np.where(ok, out[k], np.nan)
        return out


def rk_props(P, T, m, *names):
    P = np.atleast_1d(np.asarray(P, float)); T = np.atleast_1d(np.asarray(T, float)); m = np.atleast_1d(np.asarray(m, float))
    o = sf.getProp(np.array([P, T, m], dtype=object), "NH3aq_RK", sfm.defpath, *names)
    return {n: np.asarray(getattr(o, n), float) for n in names}


# ---------------------------------------------------------------------------
# Sanity of V1: pure water and handbook densities
# ---------------------------------------------------------------------------
Pc, Tc_ = np.array([0.1, 100.0, 400.0, 700.0]), np.array([240.0, 260.0, 280.0, 300.0])
a = v1_props(Pc, Tc_, [0.0])
w = sf.getProp(np.array([Pc, Tc_], dtype=object), "water1")
print("V1 at m = 0 against water1: max relative difference "
      f"rho {np.max(np.abs(a['rho'][:, :, 0] / w.rho - 1)):.1e}, sound speed {np.max(np.abs(a['vel'][:, :, 0] / w.vel - 1)):.1e}, "
      f"Cp {np.max(np.abs(a['Cp'][:, :, 0] / w.Cp - 1)):.1e}")
h = v1_props([0.1], [293.15], wt2m([10, 20]))
print(f"V1 density at 0.1 MPa, 20 C: 10 wt% {h['rho'][0, 0, 0]:.1f}, 20 wt% {h['rho'][0, 0, 1]:.1f} kg/m3")

# ---------------------------------------------------------------------------
# Figure: sound speed, density, Cp
# ---------------------------------------------------------------------------
T0, P0 = 260.0, 300.0
w_axis, P_axis, T_axis = np.linspace(0, 23, 47), np.linspace(1, 600, 61), np.linspace(240, 300, 31)
P_set, w_set = [1.0, 150.0, 300.0, 450.0], [0.0, 5.0, 10.0, 20.0]
cP = plt.get_cmap("plasma")(np.linspace(0, 0.75, len(P_set)))
cW = plt.get_cmap("viridis_r")(np.linspace(0.1, 0.95, len(w_set)))
rows = [("vel", "Sound speed (m/s)"), ("rho", "Density (kg/m$^3$)"), ("Cp", "Cp (J/kg/K)")]
fig, axes = plt.subplots(3, 3, figsize=(14.5, 11), sharex="col")
for j in range(3):
    if j == 0:
        rk, v1 = rk_props(P_set, [T0], wt2m(w_axis), "vel", "rho", "Cp"), v1_props(P_set, [T0], wt2m(w_axis))
        xs, fam, cols, lab = w_axis, P_set, cP, lambda v: f"{v:g} MPa"
        pick = lambda arr, k: arr[k, 0, :]
    elif j == 1:
        rk, v1 = rk_props(P_axis, [T0], wt2m(w_set), "vel", "rho", "Cp"), v1_props(P_axis, [T0], wt2m(w_set))
        xs, fam, cols, lab = P_axis, w_set, cW, lambda v: f"{v:g} wt%"
        pick = lambda arr, k: arr[:, 0, k]
    else:
        rk, v1 = rk_props([P0], T_axis, wt2m(w_set), "vel", "rho", "Cp"), v1_props([P0], T_axis, wt2m(w_set))
        xs, fam, cols, lab = T_axis, w_set, cW, lambda v: f"{v:g} wt%"
        pick = lambda arr, k: arr[0, :, k]
    for i, (key, ylabel) in enumerate(rows):
        ax = axes[i, j]
        for k, (val, c) in enumerate(zip(fam, cols)):
            ax.plot(xs, pick(rk[key], k), color=c, lw=1.8, label=lab(val))
            ax.plot(xs, pick(v1[key], k), color=c, lw=1.6, ls="--")
        ax.grid(True, ls=":", c="lightgray")
        if j == 0:
            ax.set_ylabel(ylabel)
    axes[0, j].set_title(["vs composition (260 K)", "vs pressure (260 K)", "vs temperature (300 MPa)"][j], fontsize=10)
    axes[0, j].legend(frameon=False, fontsize=8, ncol=2)
    axes[2, j].set_xlabel(["NH3 (wt%)", "Pressure (MPa)", "Temperature (K)"][j])
fig.suptitle("Solid: SeaFreeze 'NH3aq_RK' (Redlich-Kister, fitted to ice liquidi; these properties NOT calibrated)\n"
             "Dashed: NH3-H2O Gibbs spline V1 already in the repository (NH3_H2O_LBF.mat), shown where it is physical", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(OUT / "compare_V1_properties.png", dpi=140)
plt.close(fig)

# where is V1 physical?
Pm, Tm_ = np.arange(0.0, 800.0, 25.0), np.arange(240.0, 301.0, 5.0)
print("\nUpper pressure (MPa) up to which V1 stays physical, by temperature (rows) and NH3 wt% (columns 2, 5, 10, 20):")
phys = v1_props(Pm, Tm_, wt2m([2, 5, 10, 20]))["physical"]
for j, T in enumerate(Tm_[::2]):
    lim = []
    for k in range(4):
        col = phys[:, 2 * j, k]
        bad = np.where(~col)[0]
        lim.append(f"{Pm[bad[0]] - 25:5.0f}" if bad.size else "  775")
    print(f"  {T:5.0f} K: " + "  ".join(lim))

# numbers
print("\nAt 260 K (NH3aq_RK / V1):")
print(f"{'P MPa':>6} {'wt%':>4} | {'vel m/s':>17} {'diff':>7} | {'rho kg/m3':>17} {'diff':>7} | {'Cp J/kg/K':>15} {'diff':>6}")
for P in (1.0, 150.0, 300.0, 450.0):
    for wt in (5.0, 10.0, 20.0):
        r_ = rk_props([P], [T0], wt2m([wt]), "vel", "rho", "Cp"); v_ = v1_props([P], [T0], wt2m([wt]))
        g = lambda d_, k: float(np.ravel(d_[k])[0])
        print(f"{P:6.0f} {wt:4.0f} | {g(r_,'vel'):8.1f} {g(v_,'vel'):8.1f} {100*(g(r_,'vel')/g(v_,'vel')-1):+6.1f}% | "
              f"{g(r_,'rho'):8.1f} {g(v_,'rho'):8.1f} {100*(g(r_,'rho')/g(v_,'rho')-1):+6.1f}% | "
              f"{g(r_,'Cp'):7.0f} {g(v_,'Cp'):7.0f} {100*(g(r_,'Cp')/g(v_,'Cp')-1):+5.1f}%")

# ---------------------------------------------------------------------------
# Water activity and ice liquidus from V1 (does V1 also describe the phase equilibria?)
# ---------------------------------------------------------------------------
print("\nWater activity at 260 K (NH3aq_RK / V1 / ideal):")
for P in (1.0, 300.0):
    for wt in (5.0, 10.0, 20.0):
        m = float(wt2m(wt))
        v = v1_props([P], [T0], [0.0, m])["muw"][0, 0]
        aw_v1 = np.exp((v[1] - v[0]) / (R * T0))
        aw_rk = float(np.ravel(rk_props([P], [T0], [m], "aw")["aw"])[0])
        print(f"  {P:5.0f} MPa {wt:4.0f} wt%: {aw_rk:.4f} / {aw_v1:.4f} / {NW / (NW + m):.4f}")


def liquidus(muw_fn, ice, P, lo, hi):
    def f(T):
        pt = np.empty((1,), dtype=object); pt[0] = (P, T)
        return sf.getProp(pt, ice, sfm.defpath, "G").G[0] * MW - muw_fn(P, T)
    Tg = np.arange(lo, hi, 2.0)
    r = np.array([f(t) for t in Tg])
    ok = np.isfinite(r[:-1]) & np.isfinite(r[1:]) & (np.signbit(r[:-1]) != np.signbit(r[1:]))
    if not ok.any():
        return np.nan
    i = np.where(ok)[0][-1]
    return brentq(f, Tg[i], Tg[i + 1])


with open(REPO / "SourceData" / "NH3aq" / "liquidus_data.csv", newline="") as fh:
    data = list(csv.DictReader(line for line in fh if not line.startswith("#")))
res = {"NH3aq_RK": [], "V1": []}
n_in = 0
for row in data:
    x = float(row["x_NH3"]); m = x / (1 - x) * NW; P = float(row["P_MPa"]); T = float(row["T_K"]); ice = row["ice_fitted"]
    if T < 240.0 or P > 795.0:               # V1 domain
        continue
    n_in += 1
    t_rk = liquidus(lambda p, t: float(np.ravel(rk_props([p], [t], [m], "muw")["muw"])[0]), ice, P, 240.0, 300.0)
    t_v1 = liquidus(lambda p, t: float(v1_props([p], [t], [m])["muw"][0, 0, 0]), ice, P, 240.0, 300.0)
    res["NH3aq_RK"].append(t_rk - T); res["V1"].append(t_v1 - T)
print(f"\nIce liquidus data inside the V1 domain (T >= 240 K, P <= 795 MPa): {n_in} of {len(data)} points")
for k, v in res.items():
    v = np.array(v); ok = np.isfinite(v)
    print(f"  {k:6s}: RMS {np.sqrt(np.mean(v[ok]**2)):.2f} K, mean {v[ok].mean():+.2f} K, max |dT| {np.abs(v[ok]).max():.2f} K"
          f"{'' if ok.all() else f'  ({(~ok).sum()} points without a liquidus)'}")
print("\nfigure:", OUT / "compare_V1_properties.png")
