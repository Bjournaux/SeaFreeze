"""
Density, heat capacity and excess volume implied by the revised RK model for
H2O-NH3, as a function of P, T and composition, compared with SeaFreeze NaCl(aq).

    python plot_properties.py        (after run_fit.py)

IMPORTANT: none of these quantities was fitted. The RK model was calibrated on
ice liquidus data only. What is plotted is

    V  = x1 V_water + x2 V_NH3(liq) + V_ex ,   V_ex = R T x1 x2 sum_k v_k (x1-x2)^k / 1000
    Cp = x1 Cp_water + x2 Cp_NH3(liq)           (the model has Cp_ex = 0)

Pure liquid NH3 comes from the CoolProp Gibbs energies already exported in
../For_Baptiste/outputs (differentiated numerically; CoolProp itself is not needed).
Above its melting pressure liquid NH3 is a metastable extrapolation of that EOS.
"""

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nh3_rk as m
import seafreeze.seafreeze as sf

OUT = HERE / "outputs"
NH3_SRC = HERE.parent / "For_Baptiste" / "outputs" / "Gibbs_liq_x0.0211.csv"
M_NACL = 0.058443

p = np.array([json.loads((OUT / "fit_parameters.json").read_text())["parameters"][n] for n in m.PARAM_NAMES])


# ----------------------------------------------------------------------------
# Pure liquid NH3: V and Cp from the exported CoolProp Gibbs energy
# ----------------------------------------------------------------------------

def nh3_spline():
    d = pd.read_csv(NH3_SRC, usecols=["Pressure_MPa", "Temperature_K", "G_NH3_pure"])
    Pg = np.linspace(1, 900, 200); Tg = np.linspace(223, 273, 150)
    d["iP"] = np.abs(d.Pressure_MPa.values[:, None] - Pg[None, :]).argmin(1)
    d["iT"] = np.abs(d.Temperature_K.values[:, None] - Tg[None, :]).argmin(1)
    G = np.full((Pg.size, Tg.size), np.nan)
    G[d.iP, d.iT] = d.G_NH3_pure
    keep = np.isfinite(G).all(axis=0)                 # rows dropped by the original export (T < ~233 K)
    return RectBivariateSpline(Pg, Tg[keep], G[:, keep], kx=3, ky=3, s=0), Tg[keep].min()


NH3, T_NH3_MIN = nh3_spline()


def nh3_props(P, T):
    """V [cm3/mol], Cp [J/mol/K] of liquid NH3 on broadcast (P, T)."""
    P, T = np.broadcast_arrays(np.asarray(P, float), np.asarray(T, float))
    V = NH3.ev(P.ravel(), T.ravel(), dx=1).reshape(P.shape)
    Cp = (-T.ravel() * NH3.ev(P.ravel(), T.ravel(), dy=2)).reshape(P.shape)
    return V, Cp


def water_props(P, T):
    P, T = np.broadcast_arrays(np.asarray(P, float), np.asarray(T, float))
    rho, Cp = m._sf_scatter(P.ravel(), T.ravel(), "water1", ("rho", "Cp"))
    return (m.M_H2O / rho * 1e6).reshape(P.shape), (Cp * m.M_H2O).reshape(P.shape)


def nh3_solution(P, T, w):
    """rho [kg/m3], Cp [J/kg/K], apparent-volume excess [cm3/mol NH3] for NH3 wt% w."""
    P, T, w = np.broadcast_arrays(np.asarray(P, float), np.asarray(T, float), np.asarray(w, float))
    x2 = m.wt_to_x(w); x1 = 1 - x2
    Vw, Cpw = water_props(P, T)
    Vn, Cpn = nh3_props(P, T)
    _, _, Vex = m.excess_properties(P, T, x2, p)
    Mmix = x1 * m.M_H2O + x2 * m.M_NH3
    rho = Mmix / ((x1 * Vw + x2 * Vn + Vex) * 1e-6)
    Cp = (x1 * Cpw + x2 * Cpn) / Mmix
    # same definition as SeaFreeze 'Vex' for NaCl: apparent molar volume minus its infinite-dilution value
    eps = 1e-9
    xs = np.maximum(x2, eps)
    Va = Vn + m.excess_properties(P, T, xs, p)[2] / xs
    Va0 = Vn + m.excess_properties(P, T, eps, p)[2] / eps
    return rho, Cp, Va - Va0


M_REF = 0.2    # mol/kg (1.2 wt%), dilute reference for NaCl. Not smaller: V_phi at very low
               # molality amplifies tiny density differences between NaClaq and its m = 0 limit
               # (with 0.01 mol/kg the reference shifts by ~17 cm3/mol above the 500 MPa stitch).


def nacl_grid(P, T, w):
    """
    rho, Cp and apparent-volume excess of SeaFreeze NaClaq on a (P, T, wt%) grid -> (nP, nT, nw).

    The excess is computed here from the densities, with the same definition as for NH3:
        V_phi = [(1 + m M_s)/rho - 1/rho_w] / m ,   excess = V_phi(m) - V_phi(M_REF)
    SeaFreeze's own 'Vex' output is not used: it diverges at m = 0 and is discontinuous
    across the 500 MPa low-P/high-P stitch of NaClaq.
    """
    P = np.atleast_1d(np.asarray(P, float)); T = np.atleast_1d(np.asarray(T, float))
    w = np.atleast_1d(np.asarray(w, float))
    mol = (w / 100) / (1 - w / 100) / M_NACL
    mm, inv = np.unique(np.concatenate([[0.0, M_REF], mol]), return_inverse=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        o = sf.getProp(np.array([P, T, mm], dtype=object), "NaClaq")
    rho = np.asarray(o.rho, float).reshape(P.size, T.size, mm.size)
    Cp = np.asarray(o.Cp, float).reshape(P.size, T.size, mm.size)
    with np.errstate(divide="ignore", invalid="ignore"):
        Vphi = ((1 + mm * M_NACL) / rho - 1 / rho[:, :, [0]]) / mm * 1e6      # cm3/mol
    Vex = Vphi - Vphi[:, :, [inv[1]]]
    Vex[:, :, mm < M_REF] = np.nan
    k = inv[2:]
    return [rho[:, :, k], Cp[:, :, k], Vex[:, :, k]]


def nacl_solution(P, T, w):
    return [float(v[0, 0, 0]) for v in nacl_grid([P], [T], [w])]


# ----------------------------------------------------------------------------
# Figure: rows = property, columns = vs composition / pressure / temperature
# ----------------------------------------------------------------------------

def main():
    T0, P0 = 260.0, 300.0
    w_axis = np.linspace(0, 23, 47)
    P_axis = np.linspace(1, 900, 61)
    T_axis = np.linspace(240, 273, 34)
    P_set = [1.0, 300.0, 600.0, 900.0]
    w_set = [0.0, 5.0, 10.0, 20.0]
    cP = plt.get_cmap("plasma")(np.linspace(0, 0.8, len(P_set)))
    cW = plt.get_cmap("viridis_r")(np.linspace(0.1, 0.95, len(w_set)))

    rows = [("Density [kg/m$^3$]", 0), ("Cp [J/kg/K]", 1),
            ("V$_\\phi$ - V$_\\phi^{dilute}$ [cm$^3$/mol solute]", 2)]
    fig, axes = plt.subplots(3, 3, figsize=(14, 11), sharex="col")

    for j, (xs, xlabel, family, colors, flabel) in enumerate([
            (w_axis, "Solute [wt%]", P_set, cP, lambda v: f"{v:g} MPa"),
            (P_axis, "Pressure [MPa]", w_set, cW, lambda v: f"{v:g} wt%"),
            (T_axis, "Temperature [K]", w_set, cW, lambda v: f"{v:g} wt%")]):
        if j == 0:
            nacl = nacl_grid(family, [T0], xs)        # (nP, 1, nw)
        elif j == 1:
            nacl = nacl_grid(xs, [T0], family)        # (nP, 1, nw)
        else:
            nacl = nacl_grid([P0], xs, family)        # (1, nT, nw)
        for k, (val, c) in enumerate(zip(family, colors)):
            if j == 0:
                a = nh3_solution(val, T0, xs);  b = [v[k, 0, :] for v in nacl]
            elif j == 1:
                a = nh3_solution(xs, T0, val);  b = [v[:, 0, k] for v in nacl]
            else:
                a = nh3_solution(P0, xs, val);  b = [v[0, :, k] for v in nacl]
            for i in range(3):
                skip = (i == 2 and j > 0 and val == 0)          # excess is zero by definition at 0 wt%
                if not skip:
                    axes[i, j].plot(xs, a[i], c=c, lw=1.8, label=flabel(val))
                    axes[i, j].plot(xs, b[i], c=c, lw=1.4, ls="--")
        axes[2, j].set_xlabel(xlabel)
        axes[0, j].set_title({0: f"vs composition (T = {T0:g} K)",
                              1: f"vs pressure (T = {T0:g} K)",
                              2: f"vs temperature (P = {P0:g} MPa)"}[j], fontsize=10)
        axes[0, j].legend(frameon=False, fontsize=8, ncol=2)
    for i, (lab, _) in enumerate(rows):
        axes[i, 0].set_ylabel(lab)
        for j in range(3):
            axes[i, j].grid(True, ls=":", c="lightgray")
    axes[2, 0].axhline(0, c="k", lw=0.6); axes[2, 1].axhline(0, c="k", lw=0.6); axes[2, 2].axhline(0, c="k", lw=0.6)
    fig.suptitle("Solid: H$_2$O-NH$_3$, revised RK model (NOT fitted to these properties)      "
                 "Dashed: SeaFreeze NaCl(aq)\n"
                 "Bottom row: apparent molar volume of the solute relative to its dilute value "
                 "(NH$_3$: infinite dilution; NaCl: 0.2 mol/kg)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(OUT / "properties_PTX.png", dpi=150); fig.savefig(OUT / "properties_PTX.pdf")
    plt.close(fig)

    # numbers for a quick sanity check
    print("NH3 pure liquid (from exported CoolProp G): V, Cp at 1 MPa/273 K:",
          [round(float(v), 2) for v in nh3_props(1.0, 273.0)])
    for w in (10.0, 20.0):
        r, c, v = nh3_solution(1.0, 273.0, w)
        r_id = nh3_ideal_density(1.0, 273.0, w)
        print(f"NH3 {w:g} wt%, 1 MPa, 273 K: rho = {float(r):.1f} kg/m3 (ideal mixing {r_id:.1f}), "
              f"Cp = {float(c):.0f} J/kg/K, Vphi-Vphi_inf = {float(v):+.2f} cm3/mol")
        r, c, v = nacl_solution(1.0, 273.0, w)
        print(f"NaCl {w:g} wt%, 1 MPa, 273 K: rho = {float(r):.1f} kg/m3, Cp = {float(c):.0f} J/kg/K, "
              f"Vphi-Vphi_inf = {float(v):+.2f} cm3/mol")
    print("Written:", OUT / "properties_PTX.png")


def nh3_ideal_density(P, T, w):
    x2 = float(m.wt_to_x(w)); x1 = 1 - x2
    Vw, _ = water_props(P, T); Vn, _ = nh3_props(P, T)
    return float((x1 * m.M_H2O + x2 * m.M_NH3) / ((x1 * Vw + x2 * Vn) * 1e-6))


if __name__ == "__main__":
    main()
