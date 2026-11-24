"""
Revised Redlich-Kister (RK) model for the ice liquidus of H2O-NH3 solutions.

One liquid model is used for every ice polymorph:

    G_ex = R T x1 x2 * sum_k L_k(P,T) (x1 - x2)^k ,  k = 0, 1, 2
    L_k(P,T) = a_k + b_k (T_REF/T - 1) + v_k (P - P_REF)/1000

    x1 = x(H2O), x2 = x(NH3); P in MPa, T in K; v_k is per GPa.

Pure phases come from the installed SeaFreeze package (no pre-exported tables):
    ices Ih, II, III, V, VI  -> seafreeze.getProp
    liquid water, T >= 240 K -> 'water1' (Bollengier et al. 2019)
    liquid water, T <  240 K -> constant-Cp extrapolation of 'water1' from 240 K

Outside the domain of a SeaFreeze phase the result is NaN (never clipped).
"""

import warnings

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import brentq

import seafreeze.seafreeze as sf

R = 8.314462618        # J/mol/K
M_H2O = 0.01801528     # kg/mol
M_NH3 = 0.01703052     # kg/mol

T_REF = 240.0          # K, reference temperature (lower limit of 'water1')
P_REF = 0.1            # MPa
P_SCALE = 1000.0       # v_k are expressed per GPa

ICES = ("Ih", "II", "III", "V", "VI")
PARAM_NAMES = ("a0", "a1", "a2", "b0", "b1", "b2", "v0", "v1", "v2")


# ----------------------------------------------------------------------------
# Composition
# ----------------------------------------------------------------------------

def wt_to_x(w_percent):
    """NH3 mass percent -> NH3 mole fraction."""
    w = np.asarray(w_percent, dtype=float) / 100.0
    n_nh3 = w / M_NH3
    n_h2o = (1.0 - w) / M_H2O
    return n_nh3 / (n_nh3 + n_h2o)


def x_to_wt(x):
    """NH3 mole fraction -> NH3 mass percent."""
    x = np.asarray(x, dtype=float)
    return 100.0 * x * M_NH3 / (x * M_NH3 + (1.0 - x) * M_H2O)


# ----------------------------------------------------------------------------
# SeaFreeze access
# ----------------------------------------------------------------------------

def _scatter(P, T):
    P = np.atleast_1d(np.asarray(P, dtype=float))
    T = np.atleast_1d(np.asarray(T, dtype=float))
    P, T = np.broadcast_arrays(P, T)
    PT = np.empty((P.size,), dtype=object)
    for i, (p, t) in enumerate(zip(P.ravel(), T.ravel())):
        PT[i] = (p, t)
    return PT, P.shape


def _sf_scatter(P, T, phase, fields):
    PT, shape = _scatter(P, T)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = sf.getProp(PT, phase)
    return [np.asarray(getattr(out, f), dtype=float).reshape(shape) for f in fields]


def _sf_grid(P, T, phase, fields):
    PT = np.array([np.asarray(P, float), np.asarray(T, float)], dtype=object)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = sf.getProp(PT, phase)
    return [np.asarray(getattr(out, f), dtype=float) for f in fields]


# ----------------------------------------------------------------------------
# Pure-water reference (the only place to change to swap in another EOS)
# ----------------------------------------------------------------------------

def _extrapolate_below_tref(T, G_ref, S_ref, Cp_ref):
    """Constant-Cp continuation of G and S below T_REF (all molar)."""
    dT = T - T_REF
    G = G_ref - S_ref * dT + Cp_ref * (dT - T * np.log(T / T_REF))
    S = S_ref + Cp_ref * np.log(T / T_REF)
    return G, S


def water_ref(P, T):
    """
    Molar Gibbs energy [J/mol] and entropy [J/mol/K] of the pure liquid-water
    reference at scattered (P [MPa], T [K]) points.
    """
    P = np.atleast_1d(np.asarray(P, dtype=float))
    T = np.atleast_1d(np.asarray(T, dtype=float))
    P, T = np.broadcast_arrays(P, T)
    Tc = np.maximum(T, T_REF)
    G, S, Cp = _sf_scatter(P, Tc, "water1", ("G", "S", "Cp"))
    G, S, Cp = G * M_H2O, S * M_H2O, Cp * M_H2O
    below = T < T_REF
    if np.any(below):
        Gx, Sx = _extrapolate_below_tref(T[below], G[below], S[below], Cp[below])
        G = G.copy(); S = S.copy()
        G[below], S[below] = Gx, Sx
    return G, S


def water_ref_grid(P, T):
    """Same as water_ref on a (P, T) grid; returns arrays of shape (nP, nT)."""
    P = np.asarray(P, float); T = np.asarray(T, float)
    G = np.full((P.size, T.size), np.nan)
    S = np.full((P.size, T.size), np.nan)
    hi = T >= T_REF
    if np.any(hi):
        g, s = _sf_grid(P, T[hi], "water1", ("G", "S"))
        G[:, hi], S[:, hi] = g * M_H2O, s * M_H2O
    if np.any(~hi):
        g, s, cp = _sf_grid(P, np.array([T_REF]), "water1", ("G", "S", "Cp"))
        g, s, cp = g * M_H2O, s * M_H2O, cp * M_H2O           # (nP, 1)
        Gx, Sx = _extrapolate_below_tref(T[~hi][None, :], g, s, cp)
        G[:, ~hi], S[:, ~hi] = Gx, Sx
    return G, S


# ----------------------------------------------------------------------------
# Fusion quantities
# ----------------------------------------------------------------------------

def fusion(P, T, ice):
    """
    dG = G_ice - G_water_ref  [J/mol]   (negative below the pure melting point)
    dS = S_water_ref - S_ice  [J/mol/K] (positive)
    at scattered points. NaN outside the ice EOS domain.
    """
    Gw, Sw = water_ref(P, T)
    Gi, Si = _sf_scatter(P, T, ice, ("G", "S"))
    return Gi * M_H2O - Gw, Sw - Si * M_H2O


class FusionTables:
    """dG_fus(P,T) for every ice on a fine grid, for the liquidus solver."""

    def __init__(self, P_max=1000.0, dP=2.5, T_min=195.0, T_max=305.0, dT=0.25):
        self.P = np.concatenate([[P_REF], np.arange(dP, P_max + dP / 2, dP)])
        self.T = np.arange(T_min, T_max + dT / 2, dT)
        Gw, _ = water_ref_grid(self.P, self.T)
        self.interp = {}
        for ice in ICES:
            (Gi,) = _sf_grid(self.P, self.T, ice, ("G",))
            self.interp[ice] = RegularGridInterpolator(
                (self.P, self.T), Gi * M_H2O - Gw, bounds_error=False, fill_value=np.nan)

    def dG(self, ice, P, T):
        P, T = np.broadcast_arrays(np.asarray(P, float), np.asarray(T, float))
        pts = np.column_stack([np.maximum(P.ravel(), P_REF), T.ravel()])
        return self.interp[ice](pts).reshape(P.shape)


# ----------------------------------------------------------------------------
# Redlich-Kister liquid
# ----------------------------------------------------------------------------

def full_params(values, names):
    """Expand a reduced parameter vector into the full 9-vector (zeros elsewhere)."""
    p = np.zeros(len(PARAM_NAMES))
    for v, n in zip(values, names):
        p[PARAM_NAMES.index(n)] = v
    return p


def rk_L(P, T, p):
    """L0, L1, L2 at (P, T) for the full 9-parameter vector p."""
    tau = T_REF / np.asarray(T, float) - 1.0
    pi = (np.asarray(P, float) - P_REF) / P_SCALE
    return [p[k] + p[3 + k] * tau + p[6 + k] * pi for k in range(3)]


def _f_k(x2):
    x1 = 1.0 - x2
    return [np.ones_like(x1), 4.0 * x1 - 1.0, 12.0 * x1**2 - 8.0 * x1 + 1.0]


def ln_gamma_w(P, T, x2, p):
    """ln of the activity coefficient of water (mole-fraction scale, Raoult)."""
    x2 = np.asarray(x2, float)
    L = rk_L(P, T, p)
    f = _f_k(x2)
    return x2**2 * sum(Lk * fk for Lk, fk in zip(L, f))


def ln_a_w(P, T, x2, p):
    x2 = np.asarray(x2, float)
    return np.log1p(-x2) + ln_gamma_w(P, T, x2, p)


def design_matrix(P, T, x2, names):
    """Columns d(ln gamma_w)/d(param) for the requested parameter names."""
    P = np.asarray(P, float); T = np.asarray(T, float); x2 = np.asarray(x2, float)
    tau = T_REF / T - 1.0
    pi = (P - P_REF) / P_SCALE
    f = _f_k(x2)
    cols = []
    for n in names:
        k = int(n[1])
        mult = {"a": 1.0, "b": tau, "v": pi}[n[0]]
        cols.append(x2**2 * f[k] * mult)
    return np.column_stack(cols)


def excess_properties(P, T, x2, p):
    """G_ex [J/mol], H_ex [J/mol], V_ex [cm3/mol] implied by the parameters."""
    x2 = np.asarray(x2, float); x1 = 1.0 - x2
    d = x1 - x2
    L = rk_L(P, T, p)
    G_ex = R * T * x1 * x2 * sum(L[k] * d**k for k in range(3))
    H_ex = R * T_REF * x1 * x2 * sum(p[3 + k] * d**k for k in range(3))
    V_ex = R * T * x1 * x2 * sum(p[6 + k] / P_SCALE * d**k for k in range(3))
    return G_ex, H_ex, V_ex


# ----------------------------------------------------------------------------
# Liquidus
# ----------------------------------------------------------------------------

def liquidus_T(P, x2, ice, p, tables, T_lo=None, T_hi=None):
    """
    Liquidus temperature of `ice` at pressure P for NH3 mole fraction x2.
    Returns the highest-temperature root, NaN if none inside the tables.
    """
    T = tables.T
    if T_lo is not None or T_hi is not None:
        T = T[(T >= (T_lo or -np.inf)) & (T <= (T_hi or np.inf))]

    def f(t):
        return tables.dG(ice, P, t) - R * t * ln_a_w(P, t, x2, p)

    r = f(T)
    ok = np.isfinite(r[:-1]) & np.isfinite(r[1:])
    cross = np.where(ok & (np.signbit(r[:-1]) != np.signbit(r[1:])))[0]
    if cross.size == 0:
        return np.nan
    i = cross[-1]
    return brentq(lambda t: float(f(t)), T[i], T[i + 1])


def liquidus_stable(P, x2, p, tables, ices=ICES):
    """Stable liquidus: highest liquidus temperature among the ices."""
    best_T, best_ice = np.nan, None
    for ice in ices:
        t = liquidus_T(P, x2, ice, p, tables)
        if np.isfinite(t) and (not np.isfinite(best_T) or t > best_T):
            best_T, best_ice = t, ice
    return best_T, best_ice


# ----------------------------------------------------------------------------
# Optional: total liquid Gibbs energy (needs CoolProp for pure NH3)
# ----------------------------------------------------------------------------

def liquid_gibbs_total(P, T, x2, p):
    """
    G_liq = x1 G_water_ref + x2 G_NH3 + RT(x1 ln x1 + x2 ln x2) + G_ex  [J/mol].

    Requires CoolProp for liquid NH3. NOTE: only the water activity for
    x2 <= ~0.24 is constrained by the fit; G_ex and G_liq are extrapolations
    of the RK form and carry the arbitrary CoolProp reference state of NH3.
    """
    from CoolProp.CoolProp import PropsSI   # raises ImportError if missing
    x1 = 1.0 - x2
    Gw, _ = water_ref(P, T)
    G_nh3 = PropsSI("Gmolar", "P", P * 1e6, "T|liquid", T, "HEOS::Ammonia")
    G_ex, _, _ = excess_properties(P, T, x2, p)
    return x1 * Gw[0] + x2 * G_nh3 + R * T * (x1 * np.log(x1) + x2 * np.log(x2)) + G_ex
