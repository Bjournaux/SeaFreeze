"""
Maps of density, heat capacity and sound speed returned by the experimental 'NH3aq_RK'
material, in the same layout as residuals_vs_V1.py (pressure vs NH3 content at four
temperatures).  None of these three properties is calibrated for NH3aq_RK.

    python "NH3 Delaroque/RK_revised/maps_NH3aq.py"      (from the repository root)
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "Python"))
import seafreeze as sf                 # noqa: E402
import seafreeze.seafreeze as sfm      # noqa: E402

warnings.simplefilter("ignore")
OUT = HERE / "outputs"
wt2m = lambda w: (np.asarray(w, float) / (100.0 - np.asarray(w, float))) / sfm.mNH3_kgmol

P = np.linspace(1.0, 700.0, 71)
T = np.array([245.0, 260.0, 275.0, 290.0])
W = np.linspace(0.0, 23.0, 47)
m = wt2m(W)

o = sf.getProp(np.array([P, T, m], dtype=object), "NH3aq_RK", sfm.defpath, "rho", "Cp", "vel")
val = {k: np.asarray(getattr(o, k), float) for k in ("rho", "Cp", "vel")}
liquid = np.stack([np.squeeze(sf.whichphase(np.array([P, T, np.array([max(mi, 1e-9)])], dtype=object), "NH3aq_RK")) == 0
                   for mi in m], axis=2)

rows = [("rho", "Density", "kg/m$^3$", "viridis", np.arange(850, 1301, 50)),
        ("Cp", "Heat capacity", "J/kg/K", "magma", np.arange(3400, 4601, 100)),
        ("vel", "Sound speed", "m/s", "cividis", np.arange(1300, 2801, 100))]

fig, axes = plt.subplots(3, len(T), figsize=(4.1 * len(T), 11), sharex=True, sharey=True)
for i, (key, name, unit, cmap, levels) in enumerate(rows):
    vmin, vmax = np.nanmin(val[key]), np.nanmax(val[key])
    for j, Tj in enumerate(T):
        ax = axes[i, j]
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
    cb = fig.colorbar(pc, ax=axes[i, :], pad=0.01, aspect=18)
    cb.set_label(f"{name} ({unit})")
fig.suptitle("SeaFreeze 'NH3aq_RK': density, heat capacity and sound speed (NOT calibrated for this material)\n"
             "green line: ice liquidus; the liquid is stable above it (higher NH3 content), metastable below",
             fontsize=11)
fig.savefig(OUT / "maps_NH3aq.png", dpi=140, bbox_inches="tight")
plt.close(fig)

print("Values at 1 / 300 / 700 MPa (columns) for 0, 10, 20 wt% (rows):")
iP = [0, np.argmin(np.abs(P - 300)), len(P) - 1]; iW = [0, np.argmin(np.abs(W - 10)), np.argmin(np.abs(W - 20))]
for key, name, unit, _, _ in rows:
    for j, Tj in enumerate(T[[1, 3]]):
        jj = [1, 3][j]
        print(f"\n{name} ({unit.replace('$', '').replace('^', '')}) at {Tj:g} K")
        for w_i in iW:
            print(f"  {W[w_i]:4.0f} wt%: " + "  ".join(f"{val[key][p_i, jj, w_i]:8.1f}" for p_i in iP))
print(f"\nranges on the grid: " + "; ".join(f"{name} {np.nanmin(val[k]):.0f}-{np.nanmax(val[k]):.0f}" for k, name, *_ in rows))
print("figure:", OUT / "maps_NH3aq.png")
