"""
Phase boundary lines for SeaFreeze.

Computes (P, T) equilibrium curves between two phases of H2O (or between an
ice phase and aqueous NaCl) by finding the zero contour of the Gibbs-energy
difference. Mirrors the Matlab SF_PhaseLines / SF_phase_range / SF_WPD stack.

Public API:
    phase_range(material, path=defpath)  -> PhaseRange
    phase_lines(matA, matB, m=None, P=None, T=None, segment='all',
                path=defpath)            -> PhaseLineResult or list[PhaseLineResult]
    wpd(ax=None, solute='none', m=None, path=defpath) -> matplotlib.figure.Figure

Baptiste Journaux 2026 (Python port of the Matlab SF_PhaseLines rewrite).
"""

import os
import warnings
from collections import namedtuple
from dataclasses import dataclass, field
from typing import Optional, Union, List

import numpy as np

from mlbspline import load as mlb_load
from lbftd import evalGibbs as eg

from .seafreeze import phases, defpath, _load_spline


# Use the same MW_H2O as the Matlab SF_PhaseLines/SF_WhichPhase code so that
# pure-phase Gibbs energies (multiplied by MW to convert J/kg -> J/mol) agree
# bit-for-bit between implementations during cross-validation.
MW_H2O = 0.018015268  # kg/mol


PhaseRange = namedtuple('PhaseRange', ['P', 'T', 'm'])
PhaseRange.__new__.__defaults__ = (None,)  # m defaults to None for 2D phases


@dataclass
class PhaseLineResult:
    matA: str
    matB: str
    P: np.ndarray            # MPa
    T: np.ndarray            # K
    stable: np.ndarray       # bool, same length as P/T
    segment: str             # 'all' | 'stable' | 'meta'
    triple_points: np.ndarray = field(
        default_factory=lambda: np.empty((0, 2)))   # Nx2 [P_MPa, T_K]
    m: Optional[float] = None  # mol/kg, only for NaClaq pairs


# ---------------------------------------------------------------------------
# phase_range
# ---------------------------------------------------------------------------
def phase_range(material: str, path: str = defpath) -> PhaseRange:
    """Return the knot-domain range of the spline for `material`.

    Parameters
    ----------
    material : str
        One of the names in `seafreeze.phases` (Ih, II, III, V, VI,
        VII_X_French, water1, water2, water_IAPWS95, NaClaq).
    path : str
        Path to the splines/ directory (default: package splines folder).

    Returns
    -------
    PhaseRange(P, T, m) namedtuple with each element a 2-tuple
    (lo, hi). For 2D phases m is None.
    """
    if material not in phases:
        raise ValueError(
            f"Unknown material {material!r}. Valid: {', '.join(phases.keys())}.")
    if material == 'NaClaq':
        # Stitched LP+HP: report the intersection domain for T/m (both sub-splines
        # must be valid in the blend zone), full P coverage LP_lo → HP_hi.
        # At LP-only pressures T ≥ LP_Tlo = 229.9 K is required; using the
        # intersection avoids spurious extrapolation artifacts at phase boundaries.
        sp_lp = _load_spline(path, 'NaClaq_LP')
        sp_hp = _load_spline(path, 'NaClaq_HP')
        P = (float(sp_lp['knots'][0][0]), float(sp_hp['knots'][0][-1]))
        T = (max(float(sp_lp['knots'][1][0]),  float(sp_hp['knots'][1][0])),
             min(float(sp_lp['knots'][1][-1]), float(sp_hp['knots'][1][-1])))
        m_lo = max(float(sp_lp['knots'][2][0]),  float(sp_hp['knots'][2][0]))
        m_hi = min(float(sp_lp['knots'][2][-1]), float(sp_hp['knots'][2][-1]))
        return PhaseRange(P=P, T=T, m=(m_lo, m_hi))
    sp = _load_spline(path, material)
    knots = sp['knots']
    P = (float(knots[0][0]), float(knots[0][-1]))
    T = (float(knots[1][0]), float(knots[1][-1]))
    m = None
    if len(knots) > 2:
        m = (float(knots[2][0]), float(knots[2][-1]))
    return PhaseRange(P=P, T=T, m=m)


# ---------------------------------------------------------------------------
# Pair lookup table — copied from SF_PhaseLines.m (PAIRS_TABLE), excluding
# all VII_X_French pairs by request.
# ---------------------------------------------------------------------------
# Triple-point coordinates [T_K, P_MPa] (literature values used by Matlab v1).
_TP_IhLiqIII  = (251.165, 207.593)
_TP_IhIIIII   = (238.237, 209.885)
_TP_IIIIIV    = (249.418, 355.504)
_TP_IIVVI     = (201.934, 670.840)
_TP_IIIVLiq   = (256.164, 350.110)
_TP_VVILiq    = (273.407, 634.400)
_TP_atm       = (273.150, 0.000611)


def _tps(*pts):
    """Build an Nx2 array of [P_MPa, T_K] from list of (T_K, P_MPa) tuples."""
    return np.array([[p[1], p[0]] for p in pts], dtype=float)


# Each entry: (matA, matB, var, lo, hi, triple_points)
_PAIRS = [
    ('Ih',  'water1', 'T', _TP_IhLiqIII[0], _TP_atm[0],
        _tps(_TP_IhLiqIII, _TP_atm)),
    ('Ih',  'II',     'T', 100.0,            _TP_IhIIIII[0],
        _tps(_TP_IhIIIII)),
    ('Ih',  'III',    'T', _TP_IhIIIII[0],   _TP_IhLiqIII[0],
        _tps(_TP_IhIIIII, _TP_IhLiqIII)),
    ('II',  'III',    'T', _TP_IhIIIII[0],   _TP_IIIIIV[0],
        _tps(_TP_IhIIIII, _TP_IIIIIV)),
    ('II',  'V',      'T', _TP_IIVVI[0],     _TP_IIIIIV[0],
        _tps(_TP_IIVVI, _TP_IIIIIV)),
    ('II',  'VI',     'T', 50.0,             _TP_IIVVI[0],
        _tps(_TP_IIVVI)),
    ('III', 'V',      'T', _TP_IIIIIV[0],    _TP_IIIVLiq[0],
        _tps(_TP_IIIIIV, _TP_IIIVLiq)),
    ('III', 'water1', 'T', _TP_IhLiqIII[0],  _TP_IIIVLiq[0],
        _tps(_TP_IhLiqIII, _TP_IIIVLiq)),
    ('V',   'water1', 'T', _TP_IIIVLiq[0],   _TP_VVILiq[0],
        _tps(_TP_IIIVLiq, _TP_VVILiq)),
    ('VI',  'water1', 'T', _TP_VVILiq[0],    1000.0,
        _tps(_TP_VVILiq)),
    ('V',   'VI',     'T', _TP_IIVVI[0],     _TP_VVILiq[0],
        _tps(_TP_IIVVI, _TP_VVILiq)),
    # II <-> water1 is entirely metastable (lo > hi triggers all-meta branch).
    ('II',  'water1', 'T', np.inf,           -np.inf,
        _tps(_TP_IhIIIII, _TP_IIIIIV)),
    # NaClaq pairs - whole curve marked stable (m-dependent triple points
    # out of scope for this rewrite, matching Matlab v1.1.x).
    ('Ih',  'NaClaq', 'T', -np.inf,           np.inf,
        _tps(_TP_IhLiqIII, _TP_atm)),
    ('III', 'NaClaq', 'T', -np.inf,           np.inf,
        _tps(_TP_IhLiqIII, _TP_IIIVLiq)),
    ('V',   'NaClaq', 'T', -np.inf,           np.inf,
        _tps(_TP_IIIVLiq, _TP_VVILiq)),
    ('VI',  'NaClaq', 'T', -np.inf,           np.inf,
        _tps(_TP_VVILiq)),
    ('II',  'NaClaq', 'T', -np.inf,           np.inf,
        _tps(_TP_IhIIIII)),
]


def _lookup_pair(matA, matB):
    """Symmetric lookup. Returns (entry, swapped). Raises ValueError if absent."""
    for entry in _PAIRS:
        a, b = entry[0], entry[1]
        if a == matA and b == matB:
            return entry, False
        if a == matB and b == matA:
            return entry, True
    raise ValueError(
        f"No equilibrium line defined for the pair ({matA}, {matB}).")


# ---------------------------------------------------------------------------
# phase_lines
# ---------------------------------------------------------------------------
def phase_lines(matA: str,
                matB: str,
                m: Optional[Union[float, np.ndarray]] = None,
                P: Optional[np.ndarray] = None,
                T: Optional[np.ndarray] = None,
                segment: str = 'all',
                path: str = defpath
                ) -> Union[PhaseLineResult, List[PhaseLineResult]]:
    """Equilibrium (P, T) curve between phases `matA` and `matB`.

    Parameters
    ----------
    matA, matB : str
        Material names from `seafreeze.phases`. Order is symmetric.
    m : float or array-like, optional
        Molality (mol/kg) - required when one of the phases is 'NaClaq'.
        If a vector, returns a list of results (one per molality value).
    P, T : array-like, optional
        Override the auto-built sampling grid. P in MPa, T in K, both 1D
        sorted ascending. If omitted the grid is the spline-domain
        intersection at adaptive step (Matlab parity).
    segment : {'all', 'stable', 'meta'}, optional
        Which portion of the contour to keep.
    path : str, optional
        Path to the splines/ directory (default: package splines folder).

    Returns
    -------
    PhaseLineResult, or a list of PhaseLineResult when `m` is a vector.
    """
    segment = segment.lower()
    if segment not in {'all', 'stable', 'meta'}:
        raise ValueError(
            f"segment must be 'all', 'stable', or 'meta' (got {segment!r}).")

    entry, swapped = _lookup_pair(matA, matB)
    matA_canon, matB_canon = entry[0], entry[1]
    var, lo, hi, tps = entry[2], entry[3], entry[4], entry[5]

    nacl_involved = 'NaClaq' in (matA_canon, matB_canon)

    # Validate molality argument
    if nacl_involved:
        if m is None:
            raise ValueError(
                f"Pair ({matA_canon}, {matB_canon}) involves NaClaq; "
                f"pass molality via m= (mol/kg).")
        m_iter = np.atleast_1d(np.asarray(m, dtype=float)).ravel()
        rng_nacl = phase_range('NaClaq', path)
        if (np.any(m_iter < rng_nacl.m[0]) or
                np.any(m_iter > rng_nacl.m[1])):
            raise ValueError(
                f"Molality value(s) outside NaClaq spline range "
                f"[{rng_nacl.m[0]:g}, {rng_nacl.m[1]:g}] mol/kg: {m_iter}.")
    else:
        if m is not None:
            raise ValueError(
                f"Pair ({matA_canon}, {matB_canon}) does not involve "
                f"NaClaq; do not pass m=.")
        m_iter = np.array([np.nan])  # single iteration with no molality

    # Build PT grid (same logic regardless of m)
    P_grid, T_grid = _build_grid(matA_canon, matB_canon, P, T, path)

    results = []
    for m_val in m_iter:
        Ga, Gb = _compute_surfaces(matA_canon, matB_canon,
                                   P_grid, T_grid, m_val, path)
        Z = Ga - Gb
        T_eq, P_eq = _extract_contour(T_grid, P_grid, Z)
        if T_eq.size == 0:
            tag = (f" at m={m_val:g} mol/kg"
                   if nacl_involved else "")
            raise RuntimeError(
                f"No zero-crossing of G({matA_canon})-G({matB_canon}) "
                f"found in the sampled (P, T) grid{tag}.")

        stable_mask = _classify_stable(P_eq, T_eq, var, lo, hi)
        if segment == 'stable':
            keep = stable_mask
        elif segment == 'meta':
            keep = ~stable_mask
        else:
            keep = np.ones_like(stable_mask, dtype=bool)

        result = PhaseLineResult(
            matA=matA_canon, matB=matB_canon,
            P=P_eq[keep], T=T_eq[keep],
            stable=stable_mask[keep],
            segment=segment,
            triple_points=tps,
            m=(float(m_val) if nacl_involved else None),
        )
        results.append(result)

    if len(results) == 1:
        return results[0]
    return results


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------
def _adaptive_dP(span_MPa: float) -> float:
    """Same step ladder as Matlab SF_PhaseLines: 1/2/5/20 MPa."""
    if span_MPa <= 500.0:
        return 1.0
    if span_MPa <= 1000.0:
        return 2.0
    if span_MPa <= 5000.0:
        return 5.0
    return 20.0


def _build_grid(matA, matB, P_user, T_user, path):
    rA = phase_range(matA, path)
    rB = phase_range(matB, path)
    P_lo = max(rA.P[0], rB.P[0], 0.1)
    P_hi = min(rA.P[1], rB.P[1])
    T_lo = max(rA.T[0], rB.T[0], 1.0)
    T_hi = min(rA.T[1], rB.T[1])
    if P_lo >= P_hi or T_lo >= T_hi:
        raise ValueError(
            f"Phases {matA} and {matB} have non-overlapping spline domains.")

    if P_user is None:
        dP = _adaptive_dP(P_hi - P_lo)
        # arange-equivalent inclusive of upper bound
        P = np.arange(P_lo, P_hi + dP / 2.0, dP)
    else:
        P = np.asarray(P_user, dtype=float).ravel()

    if T_user is None:
        T = np.arange(T_lo, T_hi + 0.25, 0.5)
    else:
        T = np.asarray(T_user, dtype=float).ravel()

    return P, T


def _g_pure(material, P, T, path):
    """G * MW_H2O on a (nP, nT) grid, in J/mol of H2O."""
    sp = _load_spline(path, material)
    PTm = np.array([P, T], dtype=object)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        out = eg.evalSolutionGibbsGrid(sp, PTm, 'G', allowExtrapolations=False)
    return np.asarray(out.G) * MW_H2O


def _muw_nacl(P, T, m_val, path):
    """muw on a (nP, nT) grid, in J/mol of H2O.

    Uses getProp('NaClaq', ...) which handles LP+HP stitching automatically.
    """
    from .seafreeze import getProp
    PTm = np.array([P, T, np.array([m_val], dtype=float)], dtype=object)
    out = getProp(PTm, 'NaClaq', path, 'muw')
    G = np.squeeze(np.asarray(out.muw))
    if G.ndim == 1:
        # safety: if either P or T has length 1 and squeeze collapsed wrong
        G = G.reshape(len(P), len(T))
    return G


def _compute_surfaces(matA, matB, P, T, m_val, path):
    """Return (Ga, Gb), each shape (nP, nT), in J/mol of H2O."""
    if matA == 'NaClaq':
        Ga = _muw_nacl(P, T, m_val, path)
    else:
        Ga = _g_pure(matA, P, T, path)
    if matB == 'NaClaq':
        Gb = _muw_nacl(P, T, m_val, path)
    else:
        Gb = _g_pure(matB, P, T, path)
    return Ga, Gb


def _extract_contour(T, P, Z):
    """Zero-level set of Z(P, T) using matplotlib's contour generator.

    Returns (T_eq, P_eq) as 1D arrays (concatenated segments). Empty arrays
    if no zero-crossing exists in the sampled grid.
    """
    # Lazy matplotlib import + Agg backend so this works in headless envs
    # without disturbing interactive sessions.
    import matplotlib
    if not _interactive_backend(matplotlib):
        matplotlib.use('Agg', force=False)
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            cs = ax.contour(T, P, Z, levels=[0])
        T_eq, P_eq = [], []
        # Matplotlib >=3.8: use allsegs[level_index]; one level here.
        for seg in cs.allsegs[0]:
            if seg.size == 0:
                continue
            T_eq.append(seg[:, 0])
            P_eq.append(seg[:, 1])
    finally:
        plt.close(fig)
    if not T_eq:
        return np.array([]), np.array([])
    return np.concatenate(T_eq), np.concatenate(P_eq)


def _interactive_backend(matplotlib_module) -> bool:
    """True if a display-capable backend is already in use."""
    try:
        be = matplotlib_module.get_backend().lower()
    except Exception:
        return False
    return any(k in be for k in ('qt', 'tk', 'gtk', 'macosx', 'wx', 'nbagg'))


def _classify_stable(P_eq, T_eq, var, lo, hi):
    """Mirror the SF_PhaseLines classify_stable helper."""
    if lo > hi:
        # Inverted bounds (e.g. Inf, -Inf) -> all metastable.
        return np.zeros_like(P_eq, dtype=bool)
    if not np.isfinite(lo) and not np.isfinite(hi):
        # Both non-finite (e.g. -Inf, +Inf) -> all stable.
        return np.ones_like(P_eq, dtype=bool)
    if var == 'T':
        return (T_eq >= lo) & (T_eq <= hi)
    return (P_eq >= lo) & (P_eq <= hi)


# ---------------------------------------------------------------------------
# wpd
# ---------------------------------------------------------------------------
# Default pairs drawn by wpd() -- the 12 stable-phase pairs of the canonical
# H2O phase diagram (no VII_X_French).
_WPD_PAIRS = [
    ('Ih',  'water1'),
    ('Ih',  'II'),
    ('Ih',  'III'),
    ('II',  'III'),
    ('II',  'V'),
    ('II',  'VI'),
    ('II',  'water1'),
    ('III', 'V'),
    ('III', 'water1'),
    ('V',   'water1'),
    ('V',   'VI'),
    ('VI',  'water1'),
]


def wpd(ax=None,
        solute: str = 'none',
        m: Optional[Union[float, np.ndarray]] = None,
        show_meta: bool = True,
        phase_labels: bool = False,
        path: str = defpath):
    """Draw the H2O Water Phase Diagram.

    Parameters
    ----------
    ax : matplotlib Axes, optional
        Target axes; if omitted a new figure is created.
    solute : {'none', 'NaCl', 'NaClaq'}
        If 'NaCl' (or 'NaClaq'), overlay Ih/III/V/VI -- NaClaq melting
        curves at the molalities given by `m`.
    m : float or array-like, optional
        Molality value(s) for the NaCl overlay. Required when solute is set.
    show_meta : bool, optional
        If True (default), plot metastable extensions as dashed gray lines.
        Set to False to show only stable phase boundaries.
    phase_labels : bool, optional
        If True, annotate each stability field with its phase name (Ih, II,
        III, V, VI, Liquid).
    path : str, optional
        Path to the splines/ directory (default: package splines folder).

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 7))
    else:
        fig = ax.figure

    # Pure-water / ice pairs
    for matA, matB in _WPD_PAIRS:
        try:
            r = phase_lines(matA, matB, segment='all', path=path)
        except (RuntimeError, ValueError):
            continue
        # Plot stable and (optionally) metastable runs separately, with NaN
        # breaks so a single plot() draws disjoint segments without bridging.
        Ps, Ts = _runs_with_breaks(r.P, r.T, r.stable)
        if Ps.size:
            ax.plot(Ps, Ts, '-', color='k', linewidth=1.4,
                    label=f'{matA}-{matB}')
        if show_meta:
            Pm, Tm = _runs_with_breaks(r.P, r.T, ~r.stable)
            if Pm.size:
                ax.plot(Pm, Tm, '--', color='gray', linewidth=1.0,
                        alpha=0.7)

    # NaClaq overlay
    if solute.lower() in ('nacl', 'naclaq'):
        if m is None:
            raise ValueError("solute set; pass molality value(s) via m=.")
        m_arr = np.atleast_1d(np.asarray(m, dtype=float)).ravel()
        cmap = plt.cm.viridis
        for ice in ('Ih', 'III', 'V', 'VI'):
            for i, mv in enumerate(m_arr):
                try:
                    r = phase_lines(ice, 'NaClaq', m=mv, path=path)
                except (RuntimeError, ValueError):
                    continue
                color = cmap(i / max(len(m_arr) - 1, 1))
                ax.plot(r.P, r.T, '-', color=color, linewidth=1.2,
                        label=f'{ice}-NaClaq m={mv:g}')

    ax.set_xlabel('Pressure (MPa)', fontsize=12)
    ax.set_ylabel('Temperature (K)', fontsize=12)
    ax.set_title('H$_2$O Phase Diagram (SeaFreeze)', fontsize=13)
    ax.set_xlim(0, 2300)
    ax.set_ylim(150, 400)
    ax.grid(True, alpha=0.3)

    # Phase-field labels
    if phase_labels:
        _label_kw = dict(fontsize=11, fontweight='bold', ha='center',
                         va='center', color='#333333',
                         bbox=dict(boxstyle='round,pad=0.2', fc='white',
                                   ec='none', alpha=0.7))
        # Positions chosen to sit at the midpoint between the two bounding
        # phase lines for each field (computed from the actual equilibrium
        # curves; see phaselines.py _PAIRS for triple-point references).
        ax.text(80,   240,  'Ih',      **_label_kw)   # low-P ice
        ax.text(330,  217,  'II',      **_label_kw)   # mid-P cold ice
        ax.text(280,  248,  'III',     **_label_kw)   # between II-III and III-melt
        ax.text(490,  247,  'V',       **_label_kw)   # between II-V and V-melt
        ax.text(1200, 255,  'VI',      **_label_kw)   # well below VI melting
        ax.text(300,  345,  'Liquid',  **_label_kw)   # low-P liquid region

    return fig


def _runs_with_breaks(P, T, mask):
    """Return concatenated P, T arrays with NaN inserted between runs of
    `mask==True`. Lets a single plot() call draw disjoint segments."""
    P = np.asarray(P)
    T = np.asarray(T)
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return np.array([]), np.array([])
    # Find runs of True
    padded = np.r_[False, mask, False]
    starts = np.flatnonzero(np.diff(padded.astype(int)) == 1)
    ends = np.flatnonzero(np.diff(padded.astype(int)) == -1)
    Ps, Ts = [], []
    for s, e in zip(starts, ends):
        Ps.append(P[s:e])
        Ts.append(T[s:e])
        Ps.append(np.array([np.nan]))
        Ts.append(np.array([np.nan]))
    # strip trailing NaN
    if Ps:
        Ps = Ps[:-1]; Ts = Ts[:-1]
    return np.concatenate(Ps), np.concatenate(Ts)
