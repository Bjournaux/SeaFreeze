"""
SF_rho2P — pressure from density and temperature.

Inverts the SeaFreeze EOS to find pressure P such that rho(P, T) == rho_target.

Public API
----------
rho2P(rho_target, T, phase, m=None, *, P0=None, tol=0.01, max_iter=30, path=defpath)

Baptiste Journaux - 2026
"""

import warnings
import numpy as np

from .seafreeze import getProp, phases, defpath
from .phaselines import phase_range


def rho2P(rho_target, T, phase, m=None, *, P0=None, tol=0.01, max_iter=30, path=defpath):
    """Invert the SeaFreeze EOS: find P such that rho(P, T) == rho_target.

    Parameters
    ----------
    rho_target : float or array-like
        Target density in kg/m³.
    T : float or array-like
        Temperature in K.  Scalar is broadcast against rho_target.
    phase : str
        Any SeaFreeze material code: 'water1', 'water2', 'water_IAPWS95',
        'Ih', 'II', 'III', 'V', 'VI', 'VII_X_French',
        'NaClaq', 'NaClaq_LP', 'NaClaq_HP', 'NaClaq_5GPa_2024'.
    m : float or array-like, optional
        Molality in mol/kg — required for NaClaq phases.
        Scalar is broadcast against rho_target.
    P0 : float or array-like, optional
        Initial pressure guess in MPa.  If omitted, one Newton step from a
        stable interior point is used.
    tol : float, optional
        Convergence tolerance in MPa (default: 0.01 MPa).
    max_iter : int, optional
        Maximum Newton-Raphson iterations before bisection fallback (default 30).
    path : str, optional
        Path to the splines directory (default: package splines).

    Returns
    -------
    numpy.ndarray
        Pressure in MPa, same shape as rho_target.  NaN where no solution was
        found (rho_target outside the phase's density range at T, or T outside
        the spline domain).

    Notes
    -----
    Method: Newton-Raphson using the isothermal bulk modulus Kt (MPa)::

        P_{n+1} = P_n + (rho_target - rho(P_n, T)) * Kt(P_n, T) / rho(P_n, T)

    Convergence is verified by a final residual check; points whose
    ``|rho_final - rho_target| > 100 * tol`` are returned as NaN.
    A bisection fallback handles stalled iterations.

    Examples
    --------
    >>> import numpy as np
    >>> import seafreeze as sf
    >>> sf.rho2P(1100.0, 300.0, 'water1')       # ≈ 300 MPa
    array(299.5...)

    >>> # Ice VI scatter
    >>> sf.rho2P([1310., 1350., 1390.], [255., 260., 265.], 'VI')

    >>> # NaClaq at 1 mol/kg
    >>> sf.rho2P(1050.0, 300.0, 'NaClaq', m=1.0)
    """
    # ---- Validate phase --------------------------------------------------------
    if phase not in phases:
        raise ValueError(
            f"Unknown phase '{phase}'. Supported: {', '.join(phases.keys())}")

    is_nacl = phase.startswith('NaClaq')

    if is_nacl and m is None:
        raise ValueError(
            f"NaClaq phase '{phase}' requires the 'm' argument (molality in mol/kg).")
    if not is_nacl and m is not None:
        warnings.warn(
            f"'m' argument provided for non-NaClaq phase '{phase}' — it will be ignored.",
            UserWarning, stacklevel=2)

    # ---- Flatten inputs to 1-D, remember original shape -----------------------
    rho_arr = np.asarray(rho_target, dtype=float)
    sz = rho_arr.shape          # original shape (may be () for scalar input)
    rho_flat = rho_arr.ravel()
    n = rho_flat.size

    T_flat = np.asarray(T, dtype=float).ravel()
    if T_flat.size == 1:
        T_flat = np.full(n, T_flat[0])
    elif T_flat.size != n:
        raise ValueError(
            f"T must be scalar or the same size as rho_target ({n}), got {T_flat.size}.")

    if is_nacl:
        m_flat = np.asarray(m, dtype=float).ravel()
        if m_flat.size == 1:
            m_flat = np.full(n, m_flat[0])
        elif m_flat.size != n:
            raise ValueError(
                f"m must be scalar or the same size as rho_target ({n}), got {m_flat.size}.")
        if np.any(m_flat < 0):
            raise ValueError("Molality must be non-negative.")
    else:
        m_flat = None

    # ---- Domain bounds ---------------------------------------------------------
    rng = phase_range(phase, path)
    P_lo, P_hi = rng.P

    # P_newton_lo: floor for Newton clipping — use P_lo itself, or a tiny
    # positive value when P_lo == 0 (EOS may diverge exactly at P = 0).
    P_newton_lo = float(P_lo) if P_lo > 0 else 1e-3

    # P_init_lo: slightly interior point used to seed the initial guess via a
    # linearised step.  Larger than P_newton_lo for EOS that are non-physical
    # near P = 0 (e.g. water2 diverges there).
    P_init_lo = max(P_newton_lo, (P_hi - P_lo) * 0.001 + P_newton_lo)

    # ---- Initial guess ---------------------------------------------------------
    if P0 is not None:
        P_n = np.asarray(P0, dtype=float).ravel()
        if P_n.size == 1:
            P_n = np.full(n, P_n[0])
        P_n = np.clip(P_n, P_newton_lo, P_hi)
    else:
        rho0, Kt0 = _eval_rho_Kt(
            np.full(n, P_init_lo), T_flat, phase, m_flat, path)
        dP0 = (rho_flat - rho0) * Kt0 / np.where(np.abs(rho0) > 0, rho0, 1.0)
        dP0 = np.where(np.isfinite(dP0), dP0, 0.0)
        P_n = np.clip(P_init_lo + dP0, P_newton_lo, P_hi)

    # ---- Newton-Raphson --------------------------------------------------------
    converged = np.zeros(n, dtype=bool)
    for _ in range(int(max_iter)):
        rho_n, Kt_n = _eval_rho_Kt(P_n, T_flat, phase, m_flat, path)
        residual = rho_flat - rho_n
        dP = residual * Kt_n / np.where(np.abs(rho_n) > 0, rho_n, 1.0)
        dP = np.where(np.isfinite(dP), dP, 0.0)
        P_n = np.clip(P_n + dP, P_newton_lo, P_hi)
        converged |= (np.abs(dP) < tol)
        if converged.all():
            break

    # ---- Bisection fallback ----------------------------------------------------
    need_bisect = ~converged & np.isfinite(rho_flat)
    if need_bisect.any():
        idx = np.where(need_bisect)[0]
        P_n[idx] = _bisect_solve(
            rho_flat[idx], T_flat[idx], phase,
            None if m_flat is None else m_flat[idx],
            P_newton_lo, P_hi, tol, 60, path)

    # ---- Final residual verification -------------------------------------------
    valid_P = np.isfinite(P_n) & (P_n >= 0)
    bad = ~valid_P | ~np.isfinite(rho_flat)   # start pessimistic

    if valid_P.any():
        idx_v = np.where(valid_P)[0]
        rho_final = _eval_rho(P_n[idx_v], T_flat[idx_v], phase,
                               None if m_flat is None else m_flat[idx_v], path)
        close = np.isfinite(rho_final) & (
            np.abs(rho_final - rho_flat[idx_v]) <= 100.0 * tol)
        bad[idx_v] = ~close

    P_out = P_n.copy()
    P_out[bad] = np.nan
    # Restore original shape (including scalar → 0-d array)
    return P_out.reshape(sz) if sz != () else P_out.reshape(sz)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _build_ptm(P_arr, T_arr, phase, m_arr):
    """Build a scatter PTm array (1-D object array of tuples) for getProp."""
    n = P_arr.size
    PTm = np.empty(n, dtype=object)
    is_nacl = phase.startswith('NaClaq')
    if is_nacl:
        for i in range(n):
            PTm[i] = (float(P_arr[i]), float(T_arr[i]), float(m_arr[i]))
    else:
        for i in range(n):
            PTm[i] = (float(P_arr[i]), float(T_arr[i]))
    return PTm


def _eval_rho_Kt(P_arr, T_arr, phase, m_arr, path):
    """Evaluate rho (kg/m³) and Kt (MPa) at scatter (P, T[, m]) points."""
    PTm = _build_ptm(P_arr, T_arr, phase, m_arr)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            out = getProp(PTm, phase, path, 'rho', 'Kt')
        rho = np.asarray(getattr(out, 'rho', np.full(P_arr.size, np.nan)), dtype=float)
        Kt  = np.asarray(getattr(out, 'Kt',  np.full(P_arr.size, np.nan)), dtype=float)
    except Exception:
        rho = np.full(P_arr.size, np.nan)
        Kt  = np.full(P_arr.size, np.nan)
    return rho, Kt


def _eval_rho(P_arr, T_arr, phase, m_arr, path):
    """Evaluate only rho (kg/m³)."""
    PTm = _build_ptm(P_arr, T_arr, phase, m_arr)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            out = getProp(PTm, phase, path, 'rho')
        return np.asarray(getattr(out, 'rho', np.full(P_arr.size, np.nan)), dtype=float)
    except Exception:
        return np.full(P_arr.size, np.nan)


def _bisect_solve(rho_tgt, T_arr, phase, m_arr, P_lo, P_hi, tol, max_iter, path):
    """Bisection on f(P) = rho(P, T) - rho_tgt = 0 for a batch of points."""
    n = rho_tgt.size
    a = np.full(n, P_lo)
    b = np.full(n, P_hi)

    fa = _eval_rho(a, T_arr, phase, m_arr, path) - rho_tgt
    fb = _eval_rho(b, T_arr, phase, m_arr, path) - rho_tgt

    no_bracket = ~np.isfinite(fa) | ~np.isfinite(fb) | (fa * fb > 0)
    if no_bracket.any():
        warnings.warn(
            f"rho2P: {no_bracket.sum()} point(s) could not be bracketed for bisection.",
            UserWarning, stacklevel=4)

    for _ in range(max_iter):
        c  = (a + b) * 0.5
        fc = _eval_rho(c, T_arr, phase, m_arr, path) - rho_tgt
        left = (fa * fc <= 0)
        b  = np.where(left,  c,  b);  fb = np.where(left,  fc, fb)
        a  = np.where(~left, c,  a);  fa = np.where(~left, fc, fa)
        if np.max(np.abs(b - a)) < tol:
            break

    result = (a + b) * 0.5
    result[no_bracket] = np.nan
    return result
