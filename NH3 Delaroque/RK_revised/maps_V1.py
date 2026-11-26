"""
Maps of density, heat capacity and sound speed from the NH3-H2O Gibbs spline that was
already in the repository (Python/mlbspline/test/NH3_H2O_LBF.mat, "V1"), in the same layout
and with the same colour scales as maps_NH3aq.py.

    python "NH3 Delaroque/RK_revised/maps_V1.py"      (from the repository root)

Grey: V1 is not physical there (it oscillates at high pressure and low temperature).
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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
wt2m = lambda w: (np.asarray(w, float) / (100.0 - np.asarray(w, float))) / sfm.mNH3_kgmol

V1 = load.loadSpline(str(REPO / "Python" / "mlbspline" / "test" / "NH3_H2O_LBF.mat"), "spG_PTm_NH3_H2O_V1")
TC = float(V1["Tc"])

P = np.linspace(1.0, 700.0, 71)
T = np.array([245.0, 260.0, 275.0, 290.0])
W = np.linspace(0.0, 23.0, 47)
m = wt2m(W)

x = np.empty(3, dtype=object); x[0], x[1], x[2] = P, np.log(T / TC), m
d = lambda *o: eval_bspline(V1, x, list(o))     # NOT mlbspline's splev-based evaluator: see spline_eval.py
Tm = T[None, :, None]
GP, GPP, Gt, Gtt, GPt = d(1, 0, 0), d(2, 0, 0), d(0, 1, 0), d(0, 2, 0), d(1, 1, 0)
GTT, GPT = (Gtt - Gt) / Tm**2, GPt / Tm                        # temperature axis is ln(T/Tc)
with np.errstate(all="ignore"):
    val = dict(rho=1e6 / GP, vel=np.sqrt(GP**2 / (GPT**2 / GTT - GPP)), Cp=-Tm * GTT)
    ok = (val["Cp"] > 2000) & (val["Cp"] < 6000) & (val["rho"] > 800) & (val["rho"] < 1400) \
        & (val["vel"] > 800) & (val["vel"] < 3500)
val = {k: np.where(ok, v, np.nan) for k, v in val.items()}

# ice liquidus of NH3aq_RK, for reference only
liquid = np.stack([np.squeeze(sf.whichphase(np.array([P, T, np.array([max(mi, 1e-9)])], dtype=object), "NH3aq_RK")) == 0
                   for mi in m], axis=2)

# same colour limits and contour levels as maps_NH3aq.py
rows = [("rho", "Density", "kg/m$^3$", "viridis", np.arange(850, 1301, 50), (848.0, 1219.0)),
        ("Cp", "Heat capacity", "J/kg/K", "magma", np.arange(3000, 4601, 100), (3377.0, 4805.0)),
        ("vel", "Sound speed", "m/s", "cividis", np.arange(1300, 2801, 100), (1204.0, 2736.0))]

fig, axes = plt.subplots(3, len(T), figsize=(4.1 * len(T), 11), sharex=True, sharey=True)
for i, (key, name, unit, cmap, levels, (vmin, vmax)) in enumerate(rows):
    for j, Tj in enumerate(T):
        ax = axes[i, j]
        ax.set_facecolor("0.8")
        Z = val[key][:, j, :].T
        pc = ax.pcolormesh(P, W, Z, cmap=cmap, vmin=vmin, vmax=vmax, shading="auto")
        cs = ax.contour(P, W, Z, levels=levels, colors="w", linewidths=0.6)
        ax.clabel(cs, fmt="%g", fontsize=6)
        ax.contour(P, W, liquid[:, j, :].T.astype(float), levels=[0.5], colors="limegreen", linewidths=1.6)
        if i == 0:
            ax.set_title(f"{Tj:g} K", fontsize=11)
        if i == 2:
            ax.set_xlabel("Pressure (MPa)")
        if j == 0:
            ax.set_ylabel(f"{name}\nNH3 (wt%)")
    cb = fig.colorbar(pc, ax=axes[i, :], pad=0.01, aspect=18, extend="both")
    cb.set_label(f"{name} ({unit})")
fig.suptitle("Earlier NH3-H2O Gibbs spline V1 (NH3_H2O_LBF.mat): density, heat capacity and sound speed\n"
             "same colour scales as the NH3aq_RK maps   |   grey: V1 not physical   |   green line: ice liquidus of NH3aq_RK, for reference",
             fontsize=11)
fig.savefig(OUT / "maps_V1.png", dpi=140, bbox_inches="tight")
plt.close(fig)

print(f"V1 physical on {100 * ok.mean():.0f} % of the grid")
iP = [0, np.argmin(np.abs(P - 300)), len(P) - 1]; iW = [0, np.argmin(np.abs(W - 10)), np.argmin(np.abs(W - 20))]
print("Values at 1 / 300 / 700 MPa (columns) for 0, 10, 20 wt% (rows):")
for key, name, unit, *_ in rows:
    for jj in (1, 3):
        print(f"\n{name} ({unit.replace('$', '').replace('^', '')}) at {T[jj]:g} K")
        for w_i in iW:
            print(f"  {W[w_i]:4.0f} wt%: " + "  ".join(f"{val[key][p_i, jj, w_i]:8.1f}" for p_i in iP))
print("\nranges where physical: " + "; ".join(f"{name} {np.nanmin(val[k]):.0f}-{np.nanmax(val[k]):.0f}" for k, name, *_ in rows))
print("figure:", OUT / "maps_V1.png")
