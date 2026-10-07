"""Cached wrappers around the SeaFreeze Python API."""

import warnings

import numpy as np
import streamlit as st
from seafreeze import getProp as _getProp, phase_range as _phase_range, phase_lines as _phase_lines
from seafreeze import helmholtz_phases as _helmholtz_phases
from seafreeze.seafreeze import defpath as _defpath
from seafreeze.phaselines import _PAIRS

# ── Monkey-patch for SeaFreeze <= 1.1.1 shear modulus bug ────────────────────
# The PyPI version's _get_shear_mod_GPa receives T as a dtype=object array
# from meshgrid, so np.sqrt(T) fails with "loop of ufunc does not support
# argument 0 of type float which has no callable sqrt method".
# Fix: ensure T is cast to float before any ufunc operations.
import seafreeze.seafreeze as _sf_mod

_orig_shear = getattr(_sf_mod, "_get_shear_mod_GPa", None)

if _orig_shear is not None:
    def _patched_shear_mod_GPa(material, T, rho_kgm3):
        T = np.asarray(T, dtype=float)
        rho_kgm3 = np.asarray(rho_kgm3, dtype=float)
        return _orig_shear(material, T, rho_kgm3)

    _sf_mod._get_shear_mod_GPa = _patched_shear_mod_GPa


def _extract_results(out):
    """Pull all numpy array attributes from a SimpleNamespace into a dict."""
    result = {}
    for attr in dir(out):
        if not attr.startswith("_"):
            val = getattr(out, attr)
            if isinstance(val, np.ndarray):
                result[attr] = val
    return result


@st.cache_data(show_spinner="Computing properties...")
def compute_properties(P, T, m, material, props, mode, rhoT=False, branch="stable"):
    """Call getProp and return a dict of {prop_name: np.ndarray}.

    Parameters
    ----------
    P, T : tuple of floats
        Pressure (MPa) — or density (kg/m^3) when rhoT — and temperature (K).
    m : tuple of floats or None
        Molality (mol/kg) values; None for non-NaClaq materials.
    material : str
        Material code (key of seafreeze.phases).
    props : tuple of str
        Property symbols to compute.  Empty tuple = compute all.
    mode : str
        'scatter' or 'grid'.
    rhoT : bool
        First coordinate is density; the pressure is returned as 'P'.
    branch : str
        Helmholtz fluids at (P, T): 'stable', 'liquid' or 'vapor'.
    """
    P_arr = np.array(P, dtype=float)
    T_arr = np.array(T, dtype=float)

    if mode == "scatter":
        if m is not None:
            m_arr = np.array(m, dtype=float)
            PTm = np.empty(len(P_arr), dtype=object)
            for i in range(len(P_arr)):
                PTm[i] = (P_arr[i], T_arr[i], m_arr[i])
        else:
            PTm = np.empty(len(P_arr), dtype=object)
            for i in range(len(P_arr)):
                PTm[i] = (P_arr[i], T_arr[i])
    else:  # grid
        # Build a 1-D length-2/3 object array of float axes by explicit
        # assignment.  Do NOT use np.array([P_arr, T_arr], dtype=object): when
        # the axes have equal length NumPy collapses that into a 2-D (2, N)
        # object array, leaving each axis as an object-dtype row.  Old/bundled
        # seafreeze then feeds those object axes into spline evaluation and
        # NumPy 2.x raises "Cannot cast array data from dtype('O') to
        # dtype('float64') ... 'safe'".
        if m is not None:
            m_arr = np.array(m, dtype=float)
            PTm = np.empty(3, dtype=object)
            PTm[0], PTm[1], PTm[2] = P_arr, T_arr, m_arr
        else:
            PTm = np.empty(2, dtype=object)
            PTm[0], PTm[1] = P_arr, T_arr

    # Always call getProp with NO selective props — let it compute all,
    # then filter the result.  This avoids a bug in getProp's selective
    # shear/Vp/Vs path where scalar vs array types mix badly on grids.
    kw = {"rhoT": True} if rhoT else {}
    if material in _helmholtz_phases and not rhoT:
        kw["branch"] = branch
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")        # out-of-range points come back as NaN
        out = _getProp(PTm, material, _defpath, **kw)
    result = _extract_results(out)

    # If specific props were requested, filter
    if props:
        result = {k: v for k, v in result.items() if k in props}

    return result


@st.cache_data
def get_phase_range(material):
    """Return (P_min, P_max), (T_min, T_max), (m_min, m_max) or None."""
    rng = _phase_range(material)
    return rng.P, rng.T, rng.m


@st.cache_data(show_spinner=False)
def get_rho_range(material):
    """(rho_min, rho_max) in kg/m^3 over the material's (P, T) range.

    Helmholtz fluids report their spline box; for a Gibbs phase the density
    is evaluated on a coarse grid over its P-T range.
    """
    rng = _phase_range(material)
    if getattr(rng, "rho", None) is not None:
        return tuple(float(v) for v in rng.rho)
    PTm = np.empty(2, dtype=object)
    PTm[0] = np.linspace(max(rng.P[0], 0.0), rng.P[1], 25)
    PTm[1] = np.linspace(max(rng.T[0], 1.0), rng.T[1], 25)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rho = np.asarray(_getProp(PTm, material, _defpath, "rho").rho, float)
    rho = rho[np.isfinite(rho) & (rho > 0)]
    return (float(rho.min()), float(rho.max())) if rho.size else (1.0, 2000.0)


def _canonical(material):
    """Map GUI material names to the names used in _PAIRS."""
    if material in ("water_Bollengier2019", "water_Brown2018", "water_IAPWS95"):
        return "water_Bollengier2019"
    if material.startswith("NaClaq"):
        return "NaClaq_Brown2026"
    return material


@st.cache_data(show_spinner="Computing phase boundaries...")
def get_stability_boundaries(material):
    """Return list of (matA, matB, P_array, T_array) for all stable boundaries of material."""
    if material in _helmholtz_phases:
        return []                 # fluid: overlaid from the full diagram (core.diagrams)
    canon = _canonical(material)
    boundaries = []
    for matA, matB, *_ in _PAIRS:
        if canon not in (matA, matB):
            continue
        if matA == "NaClaq_Brown2026" or matB == "NaClaq_Brown2026":
            continue
        if matA == "water_Bollengier2019" or matB == "water_Bollengier2019":
            pass
        try:
            res = _phase_lines(matA, matB, segment="stable")
            if res.P is not None and len(res.P) > 0:
                other = matB if matA == canon else matA
                boundaries.append((canon, other, res.P, res.T))
        except Exception:
            pass
    return boundaries


@st.cache_data(show_spinner="Computing phase line...")
def get_phase_line(matA, matB, segment="stable", m=None):
    """Return (P_array, T_array) for the equilibrium line between matA and matB."""
    try:
        res = _phase_lines(matA, matB, segment=segment, m=m)
        if res.P is not None and len(res.P) > 0:
            return res.P, res.T
    except Exception:
        pass
    return None


@st.cache_data(show_spinner="Computing phase line...")
def get_phase_line_full(matA, matB, segment="all", m=None):
    """Return (P, T, stable_mask, triple_points) for the equilibrium line.

    stable_mask is a bool array aligned with P/T (True = thermodynamically
    stable portion). triple_points is an Nx2 array of [P_MPa, T_K]. Returns
    None if no contour was found.
    """
    try:
        res = _phase_lines(matA, matB, segment=segment, m=m)
        if res.P is not None and len(res.P) > 0:
            return res.P, res.T, res.stable, res.triple_points
    except Exception:
        pass
    return None
