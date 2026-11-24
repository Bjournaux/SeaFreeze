"""
Residuals of the experimental 'NH3aq' material against the NH3-H2O Gibbs spline that was
already in the repository (Python/mlbspline/test/NH3_H2O_LBF.mat, "V1"):

    residual = (NH3aq - V1) / V1     for density, heat capacity and sound speed

    python "NH3 Delaroque/RK_revised/residuals_vs_V1.py"      (from the repository root)

V1 is shown only where it is physical (it oscillates at high pressure and low temperature).
None of these three properties is calibrated for NH3aq.
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

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
M_NH3 = sfm.mNH3_kgmol
wt2m = lambda w: (np.asarray(w, float) / (100.0 - np.asarray(w, float))) / M_NH3

V1 = load.loadSpline(str(REPO / "Python" / "mlbspline" / "test" / "NH3_H2O_LBF.mat"), "spG_PTm_NH3_H2O_V1")
TC = float(V1["Tc"])


def v1_props(P, T, m):
    x = np.empty(3, dtype=object); x[0], x[1], x[2] = P, np.log(T / TC), m
    d = lambda *o: eval_bspline(V1, x, list(o))     # NOT mlbspline's splev-based evaluator: see spline_eval.py
    Tm = T[None, :, None]
    GP, GPP, Gt, Gtt, GPt = d(1, 0, 0), d(2, 0, 0), d(0, 1, 0), d(0, 2, 0), d(1, 1, 0)
    GTT, GPT = (Gtt - Gt) / Tm**2, GPt / Tm                    # temperature axis is ln(T/Tc)
    with np.errstate(all="ignore"):
        out = dict(rho=1e6 / GP, vel=np.sqrt(GP**2 / (GPT**2 / GTT - GPP)), Cp=-Tm * GTT)
        ok = (out["Cp"] > 2000) & (out["Cp"] < 6000) & (out["rho"] > 800) & (out["rho"] < 1400) \
            & (out["vel"] > 800) & (out["vel"] < 3500)
    return out, ok


P = np.linspace(1.0, 700.0, 71)
T = np.array([245.0, 260.0, 275.0, 290.0])
W = np.linspace(0.5, 23.0, 46)
m = wt2m(W)

o = sf.getProp(np.array([P, T, m], dtype=object), "NH3aq", sfm.defpath, "rho", "Cp", "vel")
rk = {k: np.asarray(getattr(o, k), float) for k in ("rho", "Cp", "vel")}
v1, ok = v1_props(P, T, m)
liquid = np.stack([np.squeeze(sf.whichphase(np.array([P, T, np.array([mi])], dtype=object), "NH3aq")) == 0
                   for mi in m], axis=2)

res = {k: np.where(ok, 100.0 * (rk[k] / v1[k] - 1.0), np.nan) for k in rk}
rows = [("rho", "Density", 8.0), ("Cp", "Heat capacity", 25.0), ("vel", "Sound speed", 25.0)]

fig, axes = plt.subplots(3, len(T), figsize=(4.1 * len(T), 11), sharex=True, sharey=True)
for i, (key, name, lim) in enumerate(rows):
    norm = TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim)
    for j, Tj in enumerate(T):
        ax = axes[i, j]
        ax.set_facecolor("0.8")                                   # V1 not physical here
        Z = res[key][:, j, :].T
        pc = ax.pcolormesh(P, W, Z, cmap="RdBu_r", norm=norm, shading="auto")
        cs = ax.contour(P, W, np.nan_to_num(Z, nan=0.0) * np.isfinite(Z), levels=[-20, -10, -5, -2, 2, 5, 10, 20],
                        colors="k", linewidths=0.5)
        ax.clabel(cs, fmt="%g%%", fontsize=6)
        ax.contour(P, W, liquid[:, j, :].T.astype(float), levels=[0.5], colors="limegreen", linewidths=1.6)
        if i == 0:
            ax.set_title(f"{Tj:g} K", fontsize=11)
        if i == 2:
            ax.set_xlabel("Pressure (MPa)")
        if j == 0:
            ax.set_ylabel(f"{name}\nNH3 (wt%)")
    cb = fig.colorbar(pc, ax=axes[i, :], pad=0.01, aspect=18, extend="both")
    cb.set_label(f"{name}: (NH3aq - V1) / V1  (%)")
fig.suptitle("Residuals of SeaFreeze 'NH3aq' against the earlier NH3-H2O Gibbs spline (V1)\n"
             "grey: V1 not physical there   |   green line: ice liquidus (liquid is stable on the high-temperature / "
             "high-concentration side)", fontsize=11)
fig.savefig(OUT / "residuals_vs_V1.png", dpi=140, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# Numbers
# ---------------------------------------------------------------------------
print(f"V1 physical on {100 * ok.mean():.0f} % of the grid (1-700 MPa, {T[0]:g}-{T[-1]:g} K, 0.5-23 wt%)\n")


def stats(mask, label):
    print(label)
    for key, name, _ in rows:
        r = res[key][mask & ok]
        a = rk[key][mask & ok] - v1[key][mask & ok]
        unit = {"rho": "kg/m3", "Cp": "J/kg/K", "vel": "m/s"}[key]
        print(f"  {name:13s}: mean {r.mean():+6.1f} %   RMS {np.sqrt((r**2).mean()):5.1f} %   "
              f"range {r.min():+6.1f} to {r.max():+6.1f} %   (mean {a.mean():+7.1f} {unit}, n = {r.size})")


stats(np.ones_like(ok), "All points where V1 is physical:")
stats(liquid, "\nOnly where the liquid is the stable phase:")
Pm, Wm = P[:, None, None], W[None, None, :]
print("\nBy pressure band (all compositions and temperatures), mean residual %:")
print(f"  {'band':>12}  {'density':>8} {'Cp':>8} {'sound':>8}")
for lo, hi in ((0, 50), (50, 150), (150, 300), (300, 450), (450, 700)):
    k = ok & np.broadcast_to((Pm >= lo) & (Pm < hi), ok.shape)
    print(f"  {lo:4d}-{hi:<4d} MPa  " + " ".join(f"{np.nanmean(res[key][k]):+8.1f}" for key, _, _ in rows))
print("\nBy composition band, mean residual %:")
for lo, hi in ((0, 5), (5, 10), (10, 15), (15, 23.1)):
    k = ok & np.broadcast_to((Wm >= lo) & (Wm < hi), ok.shape)
    print(f"  {lo:4.0f}-{hi:<4.0f} wt%  " + " ".join(f"{np.nanmean(res[key][k]):+8.1f}" for key, _, _ in rows))
print("\nfigure:", OUT / "residuals_vs_V1.png")
