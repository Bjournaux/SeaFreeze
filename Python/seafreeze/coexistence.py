"""
Two-phase coexistence at fixed temperature for fluids that cover the vapour.

A Helmholtz fluid (water_Brown2026) spans liquid and vapour, so SeaFreeze can compute

  * the saturation (vapour-pressure) curve  G_liquid(P,T) = G_vapour(P,T)
  * the sublimation curve of an ice         G_ice(P,T)    = G_vapour(P,T)

Both are solved per temperature by Newton iteration in ln P:

    dG/dP = 1/rho   =>   d(G_A - G_B)/d ln P = P (1/rho_A - 1/rho_B)

G is per kg of H2O for every phase (ice Gibbs splines and the Helmholtz fluid
share the IAPWS-95 reference state: U = S = 0 for liquid at the triple point).

Public API
----------
saturation(T, fluid='water_Brown2026', path=defpath)          -> Coexistence(P, T, rho_A, rho_B)
sublimation(T, ice='Ih', fluid='water_Brown2026', dilute_extension=True, path=defpath)
                                                     -> Coexistence(P, T, rho_A, rho_B)

Baptiste Journaux - 2026
"""
import warnings
from collections import namedtuple

import numpy as np

from .seafreeze import getProp, defpath, helmholtz_phases, phases, _load_spline, canonical_material
from lbftd import evalHelmholtz as eh

Coexistence = namedtuple('Coexistence', ['P', 'T', 'rho_A', 'rho_B'])
Coexistence.__doc__ = """Coexistence state per temperature: P (MPa), T (K), rho_A, rho_B (kg/m^3).
saturation: A = liquid, B = vapour.  sublimation: A = ice, B = vapour."""

# Starting guesses only (the result is the model's own coexistence):
# IAPWS-95 auxiliary vapour-pressure equation (Wagner & Pruss 2002, eq. 2.5)
_TC, _PC = 647.096, 22.064
_A_SAT = (-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502)
# IAPWS R14-08 sublimation pressure of ice Ih (Wagner et al. 2011)
_TT, _PT = 273.16, 611.657e-6
_A_SUB = (-0.212144006e2, 0.273203819e2, -0.610598130e1)
_B_SUB = (0.333333333e-2, 0.120666667e1, 0.170333333e1)


def psat_iapws_aux(T):
    """IAPWS-95 auxiliary saturation pressure (MPa), Wagner & Pruss 2002 eq. 2.5."""
    T = np.asarray(T, float)
    th = 1 - T / _TC
    a = _A_SAT
    s = a[0] * th + a[1] * th ** 1.5 + a[2] * th ** 3 + a[3] * th ** 3.5 + a[4] * th ** 4 + a[5] * th ** 7.5
    return _PC * np.exp(_TC / T * s)


def psub_iapws(T):
    """IAPWS R14-08 sublimation pressure of ice Ih (MPa), 50-273.16 K."""
    T = np.asarray(T, float)
    th = T / _TT
    s = sum(a * th ** b for a, b in zip(_A_SUB, _B_SUB))
    return _PT * np.exp(s / th)


def _scatter(P, T):
    out = np.empty(P.size, dtype=object)
    for i in range(P.size):
        out[i] = (float(P[i]), float(T[i]))
    return out


def _G_rho(phase, P, T, path, branch=None):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        if branch is None:
            o = getProp(_scatter(P, T), phase, path, 'G', 'rho')
        else:
            o = getProp(_scatter(P, T), phase, path, 'G', 'rho', branch=branch)
    return np.asarray(o.G, float).ravel(), np.asarray(o.rho, float).ravel()


def _solve(fA, fB, T, P0, tol=1e-12, maxit=60):
    """Newton in ln P for G_A = G_B at each T.  fA/fB(P, T) -> (G, rho)."""
    T = np.asarray(T, float).ravel()
    lnP = np.log(np.asarray(P0, float).ravel())
    live = np.isfinite(lnP)
    rA = np.full(T.size, np.nan)
    rB = np.full(T.size, np.nan)
    for _ in range(maxit):
        if not live.any():
            break
        i = np.flatnonzero(live)
        P = np.exp(lnP[i])
        GA, rhoA = fA(P, T[i])
        GB, rhoB = fB(P, T[i])
        rA[i], rB[i] = rhoA, rhoB
        slope = P * 1e6 * (1 / rhoA - 1 / rhoB)             # d(GA-GB)/d lnP, J/kg
        with np.errstate(divide='ignore', invalid='ignore'):
            step = (GA - GB) / slope
        bad = ~np.isfinite(step)
        step = np.clip(np.where(bad, 0.0, step), -1.0, 1.0)
        lnP[i] = lnP[i] - step
        lnP[i[bad]] = np.nan
        live[i] = np.isfinite(lnP[i]) & (np.abs(step) > tol)
    P = np.exp(lnP)
    # final densities at the converged pressure
    ok = np.isfinite(P)
    if ok.any():
        _, rA[ok] = fA(P[ok], T[ok])
        _, rB[ok] = fB(P[ok], T[ok])
    return P, rA, rB


def saturation(T, fluid='water_Brown2026', path=defpath):
    """Liquid-vapour coexistence (vapour pressure) of a Helmholtz fluid.

    :param T:     temperature(s) in K, below the critical temperature
    :param fluid: Helmholtz material (default 'water_Brown2026')
    :return:      Coexistence(P [MPa], T [K], rho_A = liquid, rho_B = vapour);
                  NaN where the two branches do not both exist (T >= Tc or
                  outside the surface).
    """
    fluid = canonical_material(fluid)
    if fluid not in helmholtz_phases:
        raise ValueError(f"saturation needs a Helmholtz fluid ({', '.join(sorted(helmholtz_phases))}).")
    T = np.atleast_1d(np.asarray(T, float))
    Tq = np.minimum(T, _TC - 1e-6)
    P0 = np.where(T < _TC, psat_iapws_aux(Tq), np.nan)
    P, rl, rv = _solve(lambda P, T_: _G_rho(fluid, P, T_, path, 'liquid'),
                       lambda P, T_: _G_rho(fluid, P, T_, path, 'vapor'), T, P0)
    same = ~(np.abs(rl / rv - 1) > 1e-6)                    # branches merged (near Tc)
    P[same] = np.nan; rl[same] = np.nan; rv[same] = np.nan
    return Coexistence(P=P, T=T, rho_A=rl, rho_B=rv)


_DILUTE_WARNED = False


def _warn_dilute(Tmin, fluid):
    """Warn once per session that the dilute-vapour extension is in use."""
    global _DILUTE_WARNED
    if not _DILUTE_WARNED:
        _DILUTE_WARNED = True
        warnings.warn(
            f"Sublimation below {Tmin:.4g} K (the lowest temperature of {fluid}) uses the "
            "dilute-vapour extension: the vapour is the surface's ideal-gas part (Z = 1). "
            "Non-ideality there changes p_sub by ~1e-5 relative. Pass dilute_extension=False "
            "for NaN instead. (Shown once per session.)", UserWarning, stacklevel=3)


def sublimation(T, ice='Ih', fluid='water_Brown2026', dilute_extension=True, path=defpath):
    """Ice-vapour coexistence (sublimation pressure).

    :param T:     temperature(s) in K.  The vapour comes from the Helmholtz
                  fluid, so T must lie in its range (water_Brown2026: T >= 230 K).
    :param ice:   ice phase (default 'Ih')
    :param fluid: Helmholtz material providing the vapour (default 'water_Brown2026')
    :param dilute_extension: (default True) below the fluid's lowest
                  temperature, treat the vapour as the surface's ideal-gas
                  part alone (Z = 1).  At sublimation pressures (< 10 Pa below
                  230 K) the neglected virial terms change p_sub by ~1e-5
                  relative (validated against the NIST measurements of Bielska
                  et al. 2013, 175-253 K).  It extrapolates the fluid surface,
                  so a UserWarning is issued the first time it is used in a
                  session.  False returns NaN below the surface.
    :return:      Coexistence(P [MPa], T [K], rho_A = ice, rho_B = vapour)
    """
    ice = canonical_material(ice)
    fluid = canonical_material(fluid)
    if fluid not in helmholtz_phases:
        raise ValueError(f"sublimation needs a Helmholtz fluid ({', '.join(sorted(helmholtz_phases))}).")
    if ice not in phases or phases[ice].shear_mod_parms is None:
        raise ValueError(f"{ice!r} is not an ice phase.")
    T = np.atleast_1d(np.asarray(T, float))
    P0 = psub_iapws(np.minimum(T, _TT))
    sp = _load_spline(path, fluid)
    Tmin = eh.domain(sp)[1][0]

    def vapour(P, T_):
        G, rho = _G_rho(fluid, P, T_, path, 'vapor')
        if dilute_extension:
            lo = T_ < Tmin
            if lo.any():
                G[lo], rho[lo] = eh.ideal_gas(sp, P[lo], T_[lo])
                _warn_dilute(Tmin, fluid)
        return G, rho
    P, ri, rv = _solve(lambda P, T_: _G_rho(ice, P, T_, path), vapour, T, P0)
    return Coexistence(P=P, T=T, rho_A=ri, rho_B=rv)
