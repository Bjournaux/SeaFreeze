"""
Gibbs energy of H2O-NH3 solutions used to build the SeaFreeze 'NH3aq' spline.

Per kg of water (n_w = 1/M_w mol of H2O and m mol of NH3):

    G_kgw = G_w(P,T) + m mu_NH3(P,T) + R T [m ln x2 + n_w ln x1] + (n_w + m) g_ex
    g_ex  = R T x1 x2 sum_k L_k(P,T) (x1 - x2)^k                    (Redlich-Kister)
    L_k   = a_k + b_k (T_ref/T - 1) + v_k (P - P_ref)/1000

and the spline stores G = G_kgw / (1 + M_NH3 m), in J per kg of solution.

What is calibrated: the RK parameters were fitted to ice liquidus data only
(x_NH3 <= 0.24, 9-971 MPa, 210-270 K; see "NH3 Delaroque/RK_revised").  That
constrains the chemical potential of water, hence melting curves, muw and aw.

What is NOT calibrated: everything that depends on mu_NH3 or on derivatives of
g_ex (density, heat capacity, sound speed, mus, apparent and partial molar
properties).  mu_NH3 does not enter the chemical potential of water at all.

Pure phases
    water : SeaFreeze 'water1' for T >= 240 K; below, constant-Cp continuation
            from 240 K (G, S and Cp continuous).
    NH3   : CoolProp liquid Gibbs energies tabulated on 1-900 MPa, 233.4-273 K
            (nh3_liquid_gibbs_CoolProp.csv), continued outside that range with
            constant Cp in T and a second-order expansion in P.
"""

import json
import warnings
from pathlib import Path

import numpy as np
from scipy.interpolate import make_lsq_spline
from scipy.special import xlogy

import seafreeze as sf
import seafreeze.seafreeze as sfm

HERE = Path(__file__).resolve().parent

R = 8.314462618            # J/mol/K
M_H2O = sfm.mH2O_kgmol     # kg/mol
M_NH3 = 17.03052e-3        # kg/mol
NW = 1.0 / M_H2O           # mol of H2O per kg of water
T_REF = 240.0              # K
P_REF = 0.1                # MPa
P_SCALE = 1000.0           # v_k are per GPa

_fit = json.loads((HERE / 'rk_parameters.json').read_text())
assert _fit['T_ref_K'] == T_REF and _fit['P_ref_MPa'] == P_REF
RK = {k: float(v) for k, v in _fit['parameters'].items()}       # a0..a2, b0..b2, v0..v2
FIT_INFO = {k: _fit[k] for k in ('rms_K', 'n_points', 'x_range', 'P_range_MPa', 'T_range_K')}


# ---------------------------------------------------------------------------
# Pure water [J/kg]
# ---------------------------------------------------------------------------
def water_G(P, T):
    """Specific Gibbs energy of the liquid-water reference on a (nP, nT) grid."""
    P = np.asarray(P, float).ravel(); T = np.asarray(T, float).ravel()
    out = np.empty((P.size, T.size))
    hi = T >= T_REF

    def w1(Tv, *props):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            o = sf.getProp(np.array([P, np.asarray(Tv, float)], dtype=object), 'water1', sfm.defpath, *props)
        return [np.asarray(getattr(o, p), float) for p in props]

    if hi.any():
        out[:, hi] = w1(T[hi], 'G')[0]
    if (~hi).any():
        g, s, cp = w1([T_REF], 'G', 'S', 'Cp')                       # (nP, 1)
        t = T[~hi][None, :]
        out[:, ~hi] = g - s * (t - T_REF) + cp * ((t - T_REF) - t * np.log(t / T_REF))
    return out


# ---------------------------------------------------------------------------
# Pure liquid NH3 [J/mol]
# ---------------------------------------------------------------------------
class _SmoothTable:
    """Least-squares tensor B-spline of a gridded table, with few knots.

    An interpolating spline of the tabulated values has noisy second derivatives
    at the table edges, and those are what the continuation outside the table is
    built from.  A least-squares fit with a handful of knots is smooth and still
    reproduces the table far better than needed (residual reported by .residual).
    """

    def __init__(self, P, T, G, P_knots, T_knots, k=5):
        tP = np.r_[[P[0]] * (k + 1), P_knots, [P[-1]] * (k + 1)]
        tT = np.r_[[T[0]] * (k + 1), T_knots, [T[-1]] * (k + 1)]
        cP = make_lsq_spline(P, G, tP, k=k, axis=0).c                 # (nCP, nT)
        c = make_lsq_spline(T, cP, tT, k=k, axis=1).c                 # axis first -> (nCT, nCP)
        self.sp = {'knots': np.array([tP, tT], dtype=object), 'coefs': np.moveaxis(c, 0, 1),
                   'order': np.array([k + 1, k + 1]), 'number': np.array([tP.size - k - 1, tT.size - k - 1])}
        self.residual = float(np.abs(self(P, T) - G).max())

    def __call__(self, P, T, dx=0, dy=0):
        from mlbspline.eval import evalMultivarSpline
        x = np.empty(2, dtype=object)
        x[0], x[1] = np.asarray(P, float).ravel(), np.asarray(T, float).ravel()
        return evalMultivarSpline(self.sp, x, [int(dx), int(dy)], allowExtrapolations=True)


def _load_nh3():
    raw = np.genfromtxt(HERE / 'nh3_liquid_gibbs_CoolProp.csv', delimiter=',', comments='#', skip_header=0)
    T = raw[0, 1:]; P = raw[1:, 0]; G = raw[1:, 1:]
    fit = _SmoothTable(P, T, G,
                       P_knots=[5, 15, 30, 60, 100, 150, 220, 300, 400, 500, 600, 700, 800],
                       T_knots=[243.0, 253.0, 263.0])
    return P, T, fit


_NH3_P, _NH3_T, _NH3 = _load_nh3()
NH3_TABLE_RESIDUAL = _NH3.residual      # J/mol, smooth fit vs tabulated values
NH3_P_RANGE = (float(_NH3_P[0]), float(_NH3_P[-1]))
NH3_T_RANGE = (float(_NH3_T[0]), float(_NH3_T[-1]))


def _nh3_inside_P(P, T):
    """mu_NH3 for P inside the table and any T (constant-Cp continuation in T). Grid (nP, nT)."""
    Tc = np.clip(T, *NH3_T_RANGE)
    g = _NH3(P, Tc)                       # (nP, nT) at clamped T
    s = -_NH3(P, Tc, dy=1)
    cp = -Tc[None, :] * _NH3(P, Tc, dy=2)
    t, tc = T[None, :], Tc[None, :]
    # exact where t == tc; G, S and Cp continuous at the table edges
    return g - s * (t - tc) + cp * ((t - tc) - t * np.log(t / tc))


def nh3_mu(P, T):
    """Chemical potential of pure liquid NH3 on a (nP, nT) grid, continued outside the table."""
    P = np.asarray(P, float).ravel(); T = np.asarray(T, float).ravel()
    Pc = np.clip(P, *NH3_P_RANGE)
    out = _nh3_inside_P(Pc, T)
    outside = P != Pc
    if outside.any():
        # second-order expansion in P from the nearest table edge (one-sided differences)
        h = 1.0
        for edge, sign in ((NH3_P_RANGE[0], +1.0), (NH3_P_RANGE[1], -1.0)):
            sel = outside & (Pc == edge)
            if not sel.any():
                continue
            g0, g1, g2 = (_nh3_inside_P(np.array([edge + sign * k * h]), T)[0] for k in (0, 1, 2))
            d1 = sign * (-3 * g0 + 4 * g1 - g2) / (2 * h)
            d2 = (g0 - 2 * g1 + g2) / h**2
            dP = (P[sel] - edge)[:, None]
            out[sel] = g0[None, :] + d1[None, :] * dP + 0.5 * d2[None, :] * dP**2
    return out


# ---------------------------------------------------------------------------
# Redlich-Kister excess term
# ---------------------------------------------------------------------------
def rk_L(P, T, rk=RK):
    tau = T_REF / T - 1.0
    pi = (P - P_REF) / P_SCALE
    return [rk[f'a{k}'] + rk[f'b{k}'] * tau + rk[f'v{k}'] * pi for k in range(3)]


def ln_gamma_w(P, T, x2, rk=RK):
    """ln of the activity coefficient of water (mole-fraction scale)."""
    x1 = 1.0 - x2
    L0, L1, L2 = rk_L(P, T, rk)
    return x2**2 * (L0 + L1 * (4 * x1 - 1) + L2 * (12 * x1**2 - 8 * x1 + 1))


def muw_exact(P, T, m, rk=RK):
    """Chemical potential of water [J/mol] in closed form, on a (nP, nT, nm) grid.
    Independent of the spline and of mu_NH3: used to verify the spline."""
    P = np.asarray(P, float).ravel(); T = np.asarray(T, float).ravel(); m = np.asarray(m, float).ravel()
    Pm, Tm, mm = P[:, None, None], T[None, :, None], m[None, None, :]
    x2 = mm / (NW + mm)
    return water_G(P, T)[:, :, None] * M_H2O + R * Tm * (np.log1p(-x2) + ln_gamma_w(Pm, Tm, x2, rk))


# ---------------------------------------------------------------------------
# Solution Gibbs energy [J per kg of solution]
# ---------------------------------------------------------------------------
def G_solution(P, T, m, rk=RK):
    """G(P, T, m) with broadcastable inputs (nP,1,1), (1,nT,1), (1,1,nm)."""
    P = np.asarray(P, float); T = np.asarray(T, float); m = np.asarray(m, float)
    Pv, Tv = P.ravel(), T.ravel()
    Gw = water_G(Pv, Tv)[:, :, None]                    # J/kg water
    mu2 = nh3_mu(Pv, Tv)[:, :, None]                    # J/mol NH3
    Pm, Tm, mm = Pv[:, None, None], Tv[None, :, None], m.reshape(1, 1, -1)
    n = NW + mm
    x2 = mm / n
    x1 = 1.0 - x2
    L0, L1, L2 = rk_L(Pm, Tm, rk)
    d = x1 - x2
    g_ex = R * Tm * x1 * x2 * (L0 + L1 * d + L2 * d**2)
    mix = R * Tm * (xlogy(mm, x2) + NW * np.log(x1))
    return (Gw + mm * mu2 + mix + n * g_ex) / (1.0 + M_NH3 * mm)
