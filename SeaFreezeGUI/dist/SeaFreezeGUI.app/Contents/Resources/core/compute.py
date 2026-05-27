"""Cached wrappers around the SeaFreeze Python API."""

import numpy as np
import streamlit as st
from seafreeze import getProp as _getProp, phase_range as _phase_range, phase_lines as _phase_lines
from seafreeze.seafreeze import defpath as _defpath
from seafreeze.phaselines import _PAIRS


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
def compute_properties(P, T, m, material, props, mode):
    """Call getProp and return a dict of {prop_name: np.ndarray}.

    Parameters
    ----------
    P, T : tuple of floats
        Pressure (MPa) and temperature (K) values.
    m : tuple of floats or None
        Molality (mol/kg) values; None for non-NaClaq materials.
    material : str
        Material code (key of seafreeze.phases).
    props : tuple of str
        Property symbols to compute.  Empty tuple = compute all.
    mode : str
        'scatter' or 'grid'.
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
        if m is not None:
            m_arr = np.array(m, dtype=float)
            PTm = np.array([P_arr, T_arr, m_arr], dtype=object)
        else:
            PTm = np.array([P_arr, T_arr], dtype=object)

    # Always call getProp with NO selective props — let it compute all,
    # then filter the result.  This avoids a bug in getProp's selective
    # shear/Vp/Vs path where scalar vs array types mix badly on grids.
    out = _getProp(PTm, material, _defpath)
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


def _canonical(material):
    """Map GUI material names to the names used in _PAIRS."""
    if material in ("water1", "water2", "water_IAPWS95"):
        return "water1"
    if material.startswith("NaClaq"):
        return "NaClaq"
    return material


@st.cache_data(show_spinner="Computing phase boundaries...")
def get_stability_boundaries(material):
    """Return list of (matA, matB, P_array, T_array) for all stable boundaries of material."""
    canon = _canonical(material)
    boundaries = []
    for matA, matB, *_ in _PAIRS:
        if canon not in (matA, matB):
            continue
        if matA == "NaClaq" or matB == "NaClaq":
            continue
        if matA == "water1" or matB == "water1":
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
