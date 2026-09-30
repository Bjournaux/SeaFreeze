"""
Full water phase diagram — vapour, liquid, supercritical fluid, and the ices —
by Gibbs-energy minimisation, in (P,T) and in (rho,T).

The fluid is a Helmholtz material (water3), which spans vapour, liquid and
supercritical states; the ices are the SeaFreeze Gibbs splines (Ih, II, III,
V, VI by default).  At every grid point the stable phase is the one with the
lowest specific Gibbs energy (all phases share the IAPWS-95 reference state).

Validity masks (both on by default, see phase_map):
  * an ice competes only where its spline is physical there (rho > 0,
    K_T > 0, 0 < Cp < 2 x 9R/M): outside its fitted range a Gibbs spline can
    otherwise win spuriously far above its melting point;
  * the fluid is not used more than 40 K below the melting curve of the stable
    solid (the psiEOS 'dq2026' model: IAPWS R14-08, Datchi et al. 2000,
    Queyroux et al. 2020, French & Hamel), above the triple-point pressure —
    the surface's own validity mask (it is not water there).
Where no phase is available (e.g. the ice VII/X field when VII is not
included) the diagram is left blank.

Public API
----------
phase_map(P, T, fluid='water3', ices=ICES, path=defpath) -> PhaseMap
phase_diagram_PT(P=(1e-8, 1e5), T=(150, 1800), ...)     -> DiagramPT   (what wpd_PT draws, as data)
phase_diagram_rhoT(rho=(1e-7, 4e3), T=(150, 1800), ...) -> DiagramRhoT (what wpd_rhoT draws, as data)
property_map(diag, *props, path=defpath)                -> PropertyMap (property of the stable phase)
wpd_PT(ax=None, P=(1e-8, 1e5), T=(150, 1800), ...)      -> matplotlib Figure
wpd_rhoT(ax=None, rho=(1e-7, 4e3), T=(150, 1800), ...)  -> matplotlib Figure
melt_T_dq2026(P)                                         -> melting T (K) of the stable solid

Baptiste Journaux - 2026
"""
import warnings
from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np

from .seafreeze import getProp, defpath, _load_spline, helmholtz_phases
from .coexistence import Coexistence, saturation, sublimation
from lbftd import evalHelmholtz as eh

ICES = ('Ih', 'II', 'III', 'V', 'VI')          # add 'VII_X_French' via ices=...
_LABEL = {'Ih': 'Ih', 'II': 'II', 'III': 'III', 'V': 'V', 'VI': 'VI', 'VII_X_French': 'VII/X'}
# fixed categorical order: fluid, Ih, II, III, V, VI, VII/X (validated palette)
_COLORS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7']
_TC, _PC, _RHOC = 647.096, 22.064, 322.0


@dataclass
class PhaseMap:
    """Result of phase_map on a (P, T) grid (arrays are nP x nT)."""
    P: np.ndarray              # MPa, (nP,)
    T: np.ndarray              # K, (nT,)
    names: List[str]           # [fluid, *ices]
    G: np.ndarray              # (nphase, nP, nT) J/kg, NaN outside each spline
    rho: np.ndarray            # (nphase, nP, nT) kg/m^3
    stable: np.ndarray         # (nP, nT) index into names, -1 where nothing is defined
    rho_stable: np.ndarray     # (nP, nT) density of the stable phase


@dataclass
class DiagramPT:
    """Full phase diagram in (P, T) as data (phase_diagram_PT): what wpd_PT draws."""
    pm: PhaseMap               # stability fields on the (log P, T) grid
    boundaries: list           # (i, j, P, T): G_i = G_j segments between stable phases i < j
    saturation: Optional[Coexistence]  # vapour-liquid curve (rho_A liquid, rho_B vapour);
                               # NaN where the fluid is not stable, None outside the T range
    critical: dict             # {'P', 'T', 'rho'} of the critical point (MPa, K, kg/m^3)
    triple_points: list        # triple_points(pm)
    labels: list               # (text, P, T, kind) field labels, kind 'fluid' | 'ice'
    missing: np.ndarray        # (nP, nT) True where no phase is available


@dataclass
class DiagramRhoT:
    """Full phase diagram in (rho, T) as data (phase_diagram_rhoT): what wpd_rhoT draws."""
    pm: PhaseMap               # the (P, T) phase map it is re-drawn from
    rho: np.ndarray            # (nrho,) kg/m^3
    T: np.ndarray              # (nT,) K, = pm.T
    field: np.ndarray          # (nrho, nT) index into pm.names, two_phase in the gaps, -1 undefined
    two_phase: int             # field value of the two-phase regions (= len(pm.names))
    coexistence: list          # (i, j, rho_i, rho_j, T): coexisting densities along each P-T boundary
    saturation: Optional[Coexistence]  # on the T grid (273.16 K <= T < Tc); NaN where the fluid
                               # is not stable, None outside the T range
    critical: dict             # {'P', 'T', 'rho'} of the critical point
    tie_lines: list            # triple_points(pm): the three coexisting densities at each triple T
    labels: list               # (text, rho, T, kind), kind 'fluid' | 'ice' | 'two-phase'
    xscale: str                # 'log' | 'linear'
    gap: Optional[np.ndarray] = None  # (nrho, nT, 2) int8: the two coexisting phases
                               # (i <= j; (0, 0) = liquid + vapour) in two-phase cells, -1 elsewhere


@dataclass
class PropertyMap:
    """Property of the stable phase on a phase-diagram grid (property_map)."""
    values: Dict[str, np.ndarray]  # prop -> (nP, nT) or (nrho, nT), NaN where undefined
    phase: np.ndarray          # same shape: index into names of the phase shown, -1 undefined
                               # (rho-T: the diagram's two_phase in the two-phase regions)
    ideal_gas: np.ndarray      # True where the fluid is its dilute ideal-gas extension
    coords: str                # 'PT' or 'rhoT'
    names: List[str]           # pm.names


def _grid(P, T):
    return np.array([np.asarray(P, float), np.asarray(T, float)], dtype=object)


# ---- melting curve of the stable solid (psiEOS 'dq2026', JMB 2026) ----------
def _p_melt_Ih(T):
    th = T / 273.16
    return 611.657e-6 * (1 + 0.119539337e7 * (1 - th ** 3) + 0.808183159e5 * (1 - th ** 25.75)
                         + 0.333826860e4 * (1 - th ** 103.75))


def _p_melt_hp(T):
    T = np.asarray(T, float); p = np.full(T.shape, np.nan)
    m = (T >= 251.165) & (T < 256.164); th = T[m] / 251.165; p[m] = 208.566 * (1 - 0.299948 * (1 - th ** 60))
    m = (T >= 256.164) & (T < 273.31); th = T[m] / 256.164; p[m] = 350.1 * (1 - 1.18721 * (1 - th ** 8))
    m = (T >= 273.31) & (T < 355.0); th = T[m] / 273.31; p[m] = 632.4 * (1 - 1.07476 * (1 - th ** 4.6))
    m = (T >= 355.0) & (T <= 715.0); th = T[m] / 355.0
    p[m] = 2216.0 * np.exp(0.173683e1 * (1 - 1 / th) - 0.544606e-1 * (1 - th ** 5) + 0.806106e-7 * (1 - th ** 22))
    t7 = 715 / 355
    P715 = 2216.0 * np.exp(0.173683e1 * (1 - 1 / t7) - 0.544606e-1 * (1 - t7 ** 5) + 0.806106e-7 * (1 - t7 ** 22))
    n = np.log(45000 / P715) / np.log(1600 / 715)
    m = T > 715.0; p[m] = P715 * (T[m] / 715) ** n
    return p


def melt_T_dq2026(P):
    """Melting temperature (K) of the stable solid at P (MPa).

    Port of psiEOS.m melt_T (lbf-thermo, model 'dq2026', JMB 2026): IAPWS
    R14-08 for ice Ih/III/V/VI to 2.17 GPa; ice VII from Datchi et al. (2000)
    and superionic VII'' from Queyroux et al. (2020) to 45 GPa; above, linear
    in ln P onto French & Hamel's superionic -> fluid line.  273.16 K below
    the triple-point pressure.
    """
    P = np.asarray(P, float); shape = P.shape; P = P.ravel(); Tm = np.full(P.shape, np.nan)
    PT_si = np.array([[45, 850 * ((45 - 14.6) / 3.44 + 1) ** (1 / 4.33)], [330.6, 6000], [627.6, 7000],
                      [1034.2, 8000], [1592.9, 10000], [3734.6, 10000], [4664.8, 12000], [9311.3, 12000]])
    Tm[P < 611.657e-6] = 273.16
    ih = (P >= 611.657e-6) & (P < 208.566)
    a = np.full(ih.sum(), 251.165); b = np.full(ih.sum(), 273.16); Pi = P[ih]
    for _ in range(60):
        m = 0.5 * (a + b); up = _p_melt_Ih(m) > Pi; a = np.where(up, m, a); b = np.where(up, b, m)
    Tm[ih] = 0.5 * (a + b)
    hp = (P >= 208.566) & (P < 2170)
    a = np.full(hp.sum(), 251.165); b = np.full(hp.sum(), 2e4); Ph = P[hp]
    for _ in range(80):
        m = 0.5 * (a + b); up = _p_melt_hp(m) < Ph; a = np.where(up, m, a); b = np.where(up, b, m)
    Tm[hp] = 0.5 * (a + b)
    vii = (P >= 2170) & (P <= 45000)
    Pg = P[vii] / 1e3
    Td = 354.8 * np.maximum((Pg - 2.17) / 1.253 + 1, 1e-12) ** (1 / 3)
    xq = (Pg - 14.6) / 3.44 + 1
    Tq = np.where(xq > 0, 850 * np.maximum(xq, 0) ** (1 / 4.33), 0)
    Tm[vii] = np.maximum(Td, Tq)
    si = P > 45000
    Tm[si] = np.interp(np.log(np.minimum(P[si] / 1e3, PT_si[-1, 0])), np.log(PT_si[:, 0]), PT_si[:, 1])
    return Tm.reshape(shape)


_CP_MAX = 2 * 9 * 8.314462618 / 0.018015268     # 2 x classical 9R/M, J/kg/K


def phase_map(P, T, fluid='water3', ices=ICES, dilute_extension=True, sanity=True,
              melt_mask=True, path=defpath):
    """Stable phase of H2O on a (P, T) grid by Gibbs-energy minimisation.

    :param P, T:   1-D arrays (MPa, K)
    :param fluid:  Helmholtz material for vapour/liquid/supercritical (water3)
    :param ices:   ice phases to include (default Ih, II, III, V, VI)
    :param dilute_extension: below the fluid's lowest temperature (230 K),
                   use the fluid's ideal-gas part at P < 1e-4 MPa so the
                   vapour field and sublimation line continue to low T
    :param sanity: ignore an ice where its spline is unphysical
                   (rho <= 0, K_T <= 0, Cp <= 0 or Cp > 2 x 9R/M)
    :param melt_mask: ignore the fluid more than 40 K below melt_T_dq2026(P)
                   at P above the triple-point pressure (psiEOS validity mask)
    """
    if fluid not in helmholtz_phases:
        raise ValueError(f'fluid must be a Helmholtz material ({", ".join(sorted(helmholtz_phases))}).')
    P = np.asarray(P, float).ravel(); T = np.asarray(T, float).ravel()
    names = [fluid] + list(ices)
    G = np.full((len(names), P.size, T.size), np.nan)
    rho = np.full_like(G, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        o = getProp(_grid(P, T), fluid, path, 'G', 'rho')
        G[0], rho[0] = o.G, o.rho
        if melt_mask:
            Tm_melt = melt_T_dq2026(P)[:, None]
            bad = (P[:, None] > 611.657e-6) & (T[None, :] < Tm_melt - 40.0)
            G[0][bad] = np.nan
        if dilute_extension:
            sp = _load_spline(path, fluid)
            Tmin = eh.domain(sp)[1][0]
            Pm, Tm = np.meshgrid(P, T, indexing='ij')
            lo = (Tm < Tmin) & (Pm < 1e-4)
            if lo.any():
                G[0][lo], rho[0][lo] = eh.ideal_gas(sp, Pm[lo], Tm[lo])
        for k, ice in enumerate(ices, start=1):
            o = getProp(_grid(P, T), ice, path, 'G', 'rho', 'Cp', 'Kt')
            g = np.asarray(o.G, float)
            g[g == 0] = np.nan                                   # out-of-range sentinel
            if sanity:
                cp, kt, r = np.asarray(o.Cp), np.asarray(o.Kt), np.asarray(o.rho)
                g[~((r > 0) & (kt > 0) & (cp > 0) & (cp < _CP_MAX))] = np.nan
            G[k], rho[k] = g, o.rho
    allnan = np.all(np.isnan(G), axis=0)
    stable = np.where(allnan, -1, np.nanargmin(np.where(allnan[None], 0.0, G), axis=0))
    rs = np.take_along_axis(rho, np.maximum(stable, 0)[None], axis=0)[0]
    rs[stable < 0] = np.nan
    return PhaseMap(P=P, T=T, names=names, G=G, rho=rho, stable=stable, rho_stable=rs)


def _boundaries(pm):
    """(i, j, P_line, T_line) segments where phases i and j coexist."""
    from matplotlib.figure import Figure        # contouring only: no pyplot, no backend
    out = []
    n = len(pm.names)
    ax = Figure().subplots()
    for i in range(n):
        for j in range(i + 1, n):
            pair = (pm.stable == i) | (pm.stable == j)
            if not ((pm.stable == i).any() and (pm.stable == j).any()):
                continue
            Z = np.where(pair, pm.G[i] - pm.G[j], np.nan)
            if np.all(np.isnan(Z)):
                continue
            cs = ax.contour(np.log10(pm.P), pm.T, Z.T, levels=[0.0])
            for seg in cs.allsegs[0]:
                if len(seg) > 2:
                    out.append((i, j, 10 ** seg[:, 0], seg[:, 1]))
    return out


def _scatter1(P, T):
    pts = np.empty(1, dtype=object)
    pts[0] = (float(P), float(T))
    return pts


def _scatter(X, T):
    """getProp scatter set: 1-D object array of (X, T) tuples."""
    pts = np.empty(X.size, dtype=object)
    for q, xt in enumerate(zip(X.tolist(), T.tolist())):
        pts[q] = xt
    return pts


def _GSrho(name, P, T, fluid, path):
    """G, S, rho of one phase at one (P, T) (fluid: stable branch)."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        o = getProp(_scatter1(P, T), name, path, 'G', 'S', 'rho')
    get = lambda k: float(np.ravel(getattr(o, k, [np.nan]))[0])   # out of range: no fields
    return get('G'), get('S'), get('rho')


def triple_points(pm, fluid='water3', path=defpath):
    """Triple points of a phase map, refined by Newton on G_a = G_b = G_c.

    Candidates are the 2x2 grid blocks of pm.stable holding three distinct
    phases; each is refined in (P, T) with dG/dP = 1/rho, dG/dT = -S.  The
    ice Ih - liquid - vapour point (the fluid counts once in the map) is
    added where the saturation and sublimation curves cross.

    :return: list of dicts {phases: (a, b, c), labels, P [MPa], T [K],
             rho: (rho_a, rho_b, rho_c) [kg/m^3]}
    """
    S = pm.stable
    blocks = np.stack([S[:-1, :-1], S[1:, :-1], S[:-1, 1:], S[1:, 1:]], -1)
    cands = {}
    srt = np.sort(blocks, axis=-1)
    ndist = 1 + (np.diff(srt, axis=-1) != 0).sum(-1)
    for i, j in zip(*np.nonzero((ndist >= 3) & (srt[..., 0] >= 0))):
        key = tuple(np.unique(srt[i, j])[:3])
        cands.setdefault(key, []).append((np.log(pm.P[i]), pm.T[j]))
    out = []
    for key, pts in cands.items():
        lnP, T = np.median(np.array(pts), axis=0)
        P = np.exp(lnP)
        names = [pm.names[k] for k in key]
        ok = False
        for _ in range(30):
            v = [_GSrho(nm, P, T, fluid, path) for nm in names]
            if not np.all(np.isfinite(v)):
                break
            (Ga, Sa, ra), (Gb, Sb, rb), (Gc, Sc, rc) = v
            r = np.array([Ga - Gb, Ga - Gc])
            J = np.array([[1e6 * (1 / ra - 1 / rb), -(Sa - Sb)],
                          [1e6 * (1 / ra - 1 / rc), -(Sa - Sc)]])
            try:
                dP, dT = np.linalg.solve(J, -r)
            except np.linalg.LinAlgError:
                break
            dP = np.clip(dP, -0.5 * P, 0.5 * P); dT = np.clip(dT, -10, 10)
            P, T = P + dP, T + dT
            if abs(dP) < 1e-9 * max(P, 1e-6) and abs(dT) < 1e-9:
                ok = True
                break
        if ok:
            rho = tuple(_GSrho(nm, P, T, fluid, path)[2] for nm in names)
            out.append(dict(phases=names, labels=[_name(nm, k) for nm, k in zip(names, key)],
                            P=P, T=T, rho=rho))
    # ice Ih - liquid - vapour: saturation meets sublimation
    if 'Ih' in pm.names and pm.T.min() < 273.16 < pm.T.max():
        f = lambda T: (np.log(saturation([T], fluid=fluid, path=path).P[0])
                       - np.log(sublimation([T], fluid=fluid, path=path).P[0]))
        a, b = 265.0, 280.0
        fa = f(a)
        for _ in range(50):
            m = 0.5 * (a + b); fm = f(m)
            if np.sign(fm) == np.sign(fa):
                a, fa = m, fm
            else:
                b = m
        Tt = 0.5 * (a + b)
        st = saturation([Tt], fluid=fluid, path=path)
        rI = _GSrho('Ih', st.P[0], Tt, fluid, path)[2]
        out.append(dict(phases=['Ih', fluid, fluid], labels=['Ih', 'L', 'V'], P=float(st.P[0]), T=Tt,
                        rho=(rI, float(st.rho_A[0]), float(st.rho_B[0]))))
    return out


def _name(nm, k):
    return 'L' if k == 0 else _LABEL.get(nm, nm)


def _tint(hexcolor, a=0.30):
    c = np.array([int(hexcolor[k:k + 2], 16) for k in (1, 3, 5)]) / 255.0
    return tuple(1 - a * (1 - c))


def _style(ax):
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)


# ---- the diagrams as data ----------------------------------------------------
def _window_sat(T, fluid, path):
    """Saturation curve on _sat_T temperatures inside the T window, or None."""
    lo, hi = max(273.16, T[0]), min(_TC, T[1])
    if lo >= hi:
        return None
    Ts = _sat_T(lo)
    Ts = Ts[Ts <= hi]
    return saturation(Ts, fluid=fluid, path=path) if Ts.size else None


def _mask_coex(c, keep):
    nan = lambda a: np.where(keep, a, np.nan)
    return Coexistence(P=nan(c.P), T=c.T, rho_A=nan(c.rho_A), rho_B=nan(c.rho_B))


def _in_window(tps, P, T):
    return [tp for tp in tps if (P is None or P[0] <= tp['P'] <= P[1]) and T[0] <= tp['T'] <= T[1]]


def phase_diagram_PT(P=(1e-8, 1e5), T=(150.0, 1800.0), nP=500, nT=420, fluid='water3',
                     ices=ICES, path=defpath):
    """Full H2O phase diagram in (P, T) as data: what wpd_PT draws.

    Stability fields by Gibbs-energy minimisation (phase_map) on a log-P grid;
    boundaries are the G_i = G_j contours between neighbouring stable phases;
    the saturation curve (stable part), the critical point, the triple points
    inside the window and the field labels.

    :param P, T:   (min, max) window, MPa and K; nP log-spaced, nT linear points
    :return:       DiagramPT
    """
    Pg = np.geomspace(P[0], P[1], nP)
    Tg = np.linspace(T[0], T[1], nT)
    pm = phase_map(Pg, Tg, fluid, ices, path=path)
    sat = _window_sat(T, fluid, path)
    if sat is not None:
        sat = _mask_coex(sat, pm_stable_fluid(pm, sat.P, sat.T))       # stable branch only
    return DiagramPT(pm=pm, boundaries=_boundaries(pm), saturation=sat,
                     critical=dict(P=_PC, T=_TC, rho=_RHOC),
                     triple_points=_in_window(triple_points(pm, fluid, path), P, T),
                     labels=_labels_PT(pm), missing=pm.stable < 0)


def _labels_PT(pm):
    """(text, P, T, kind) at the median of each field (fluid split by kind)."""
    lPm, Tm = np.meshgrid(np.log10(pm.P), pm.T, indexing='ij')
    N = pm.stable.size
    out = []
    kind = _fluid_kind(pm)
    for lab, k in (('vapour', 0), ('liquid', 1), ('supercritical fluid', 2)):
        m = kind == k
        if m.sum() > max(20, 1e-3 * N):
            out.append((lab, 10 ** np.median(lPm[m]), np.median(Tm[m]), 'fluid'))
    for idx, name in enumerate(pm.names[1:], start=1):
        m = pm.stable == idx
        if m.sum() > max(4, 1.5e-4 * N):
            out.append((_LABEL.get(name, name), 10 ** np.median(lPm[m]), np.median(Tm[m]), 'ice'))
    return out


def phase_diagram_rhoT(rho=(1e-7, 4e3), T=(150.0, 1800.0), P=(1e-10, 1e5), nP=1200, nT=420,
                       nrho=800, xscale='log', fluid='water3', ices=ICES, path=defpath):
    """Full H2O phase diagram in (rho, T) as data: what wpd_rhoT draws.

    The (P, T) phase map is re-drawn in density: along every isotherm the
    density of the stable phase increases with P; where it jumps (between two
    phases, or across the vapour-liquid saturation) the density gap is a
    two-phase region -- the vapour-liquid dome, sublimation (ice + V), melting
    (ice + L) and ice-ice regions -- bounded by the coexisting densities.

    :param rho, T: (min, max) window, kg/m^3 and K
    :param P:      pressure span of the underlying (P, T) map (nP log-spaced points)
    :param xscale: 'log' (density grid log-spaced, shows the vapour) or 'linear'
    :return:       DiagramRhoT
    """
    Pg = np.geomspace(P[0], P[1], nP)
    Tg = np.linspace(T[0], T[1], nT)
    pm = phase_map(Pg, Tg, fluid, ices, path=path)
    rq = np.geomspace(rho[0], rho[1], nrho) if xscale == 'log' else np.linspace(rho[0], rho[1], nrho)
    n = len(pm.names)
    TWO = n
    # saturation densities bound the vapour-liquid dome exactly
    Tsat = Tg[(Tg >= 273.16) & (Tg < _TC)]
    sat = None
    if Tsat.size:
        # solve on ~100 temperatures clustered toward Tc, interpolate onto the grid
        Ts = _sat_T(273.16)
        s0 = saturation(Ts, fluid=fluid, path=path)
        ok = np.isfinite(s0.P)
        sat = Coexistence(P=np.exp(np.interp(Tsat, Ts[ok], np.log(s0.P[ok]))), T=Tsat,
                          rho_A=np.interp(Tsat, Ts[ok], s0.rho_A[ok]),
                          rho_B=np.exp(np.interp(Tsat, Ts[ok], np.log(s0.rho_B[ok]))))
    img = np.full((nrho, nT), -1)
    gap = np.full((nrho, nT, 2), -1, dtype=np.int8)
    pairs = {}
    for j in range(nT):
        st = pm.stable[:, j]; rs = pm.rho_stable[:, j]
        ok = (st >= 0) & np.isfinite(rs)
        if ok.sum() < 2:
            continue
        st, rs = st[ok], np.maximum.accumulate(rs[ok])
        k = np.searchsorted(rs, rq)
        inside = (k > 0) & (k < rs.size)
        kk = np.clip(k, 1, rs.size - 1)
        a, b = st[kk - 1], st[kk]
        single = inside & (a == b)
        img[single, j] = a[single]
        dual = inside & (a != b)
        img[dual, j] = TWO
        gap[dual, j, 0] = np.minimum(a, b)[dual]
        gap[dual, j, 1] = np.maximum(a, b)[dual]
        for q in np.flatnonzero(dual):
            key = (int(min(a[q], b[q])), int(max(a[q], b[q])))
            pairs.setdefault(key, []).append((rq[q], Tg[j]))
        # vapour-liquid dome from the saturation densities
        if sat is not None and Tg[j] in Tsat:
            i = int(np.flatnonzero(Tsat == Tg[j])[0])
            rv, rl = sat.rho_B[i], sat.rho_A[i]
            if np.isfinite(rv) and np.isfinite(rl):
                dome = (rq > rv) & (rq < rl) & (img[:, j] == 0)
                img[dome, j] = TWO
                gap[dome, j] = 0
                for q in np.flatnonzero(dome):
                    pairs.setdefault((0, 0), []).append((rq[q], Tg[j]))
    # coexisting densities along every P-T boundary; the ice - fluid boundary
    # passes through the triple point, where the fluid density jumps from the
    # vapour to the liquid: break the curves there
    coex = [(i, j, *_split_jumps(_rho_on(pm, i, Pl, Tl, fluid, path),
                                 _rho_on(pm, j, Pl, Tl, fluid, path), Tl))
            for i, j, Pl, Tl in _boundaries(pm)]
    if sat is not None:
        sat = _mask_coex(sat, pm_stable_fluid(pm, sat.P, Tsat))
    return DiagramRhoT(pm=pm, rho=rq, T=Tg, field=img, two_phase=TWO, coexistence=coex,
                       saturation=sat, critical=dict(P=_PC, T=_TC, rho=_RHOC),
                       tie_lines=_in_window(triple_points(pm, fluid, path), None, T),
                       labels=_labels_rhoT(pm, img, pairs, rq, Tg, xscale), xscale=xscale, gap=gap)


def _split_jumps(ri, rj, T, factor=20.0):
    """Insert NaN between consecutive points where either density changes by > factor."""
    ri, rj, T = (np.asarray(a, float) for a in (ri, rj, T))
    with np.errstate(divide='ignore', invalid='ignore'):
        jump = np.zeros(max(T.size - 1, 0), bool)
        for r in (ri, rj):
            q = np.abs(np.log(r[1:] / r[:-1]))
            jump |= q > np.log(factor)
    k = np.flatnonzero(jump) + 1
    if k.size == 0:
        return ri, rj, T
    return (np.insert(ri, k, np.nan), np.insert(rj, k, np.nan), np.insert(T, k, np.nan))


def rhoT_labels(d, xscale=None):
    """Field and two-phase labels of a DiagramRhoT, placed for a 'log' or
    'linear' density axis (default: the diagram's own grid spacing).

    The two-phase regions are labelled from d.gap (the coexisting phases of
    each two-phase cell).
    """
    xscale = xscale or d.xscale
    pairs = {}
    if d.gap is not None:
        R, TT = np.meshgrid(d.rho, d.T, indexing='ij')
        g = d.gap.astype(int)
        code = np.where(g[..., 0] >= 0, g[..., 0] * 64 + g[..., 1], -1)
        for c in np.unique(code[code >= 0]):
            m = code == c
            pairs[(int(c // 64), int(c % 64))] = list(zip(R[m], TT[m]))
    return _labels_rhoT(d.pm, d.field, pairs, d.rho, d.T, xscale)


def _wmedian(v, w):
    """Weighted median."""
    o = np.argsort(v)
    c = np.cumsum(w[o])
    return v[o][np.searchsorted(c, 0.5 * c[-1])]


def _labels_rhoT(pm, img, pairs, rq, Tg, xscale):
    """(text, rho, T, kind) for the fields and the two-phase regions.

    Positions are area-weighted medians in the displayed coordinates
    (xscale), whatever the spacing of the density grid.
    """
    fx = np.log10 if xscale == 'log' else (lambda v: v)
    ix = (lambda v: 10 ** v) if xscale == 'log' else (lambda v: v)
    xw = np.abs(np.gradient(fx(rq)))                      # displayed width of each density cell
    X, TT = np.meshgrid(fx(rq), Tg, indexing='ij')
    W = np.broadcast_to(xw[:, None], X.shape)
    N = img.size
    area = lambda m: W[m].sum() / xw.sum() * Tg.size      # in "columns x rows" units
    out = []
    for idx, name in enumerate(pm.names):
        m = img == idx
        if m.sum() < max(4, 1.2e-4 * N):
            continue
        if idx == 0:
            for lab, sel in (('vapour', (TT < _TC) & (X < fx(_RHOC))),
                             ('liquid', (TT < _TC) & (X >= fx(_RHOC))),
                             ('supercritical fluid', TT >= _TC + 50)):
                mm = m & sel
                if mm.sum() > max(15, 4.5e-4 * N):
                    out.append((lab, ix(_wmedian(X[mm], W[mm])), _wmedian(TT[mm], W[mm]), 'fluid'))
        else:
            out.append((_LABEL.get(name, name), ix(_wmedian(X[m], W[m])), _wmedian(TT[m], W[m]), 'ice'))
    wmap = dict(zip(rq, xw))
    for (i, jph), pts in pairs.items():
        pts = np.array(pts)
        if len(pts) < max(6, 1.8e-4 * N):
            continue
        if i == 0 and jph == 0:
            lab = 'L + V'
        elif i == 0:
            fl = 'V' if np.median(pts[:, 0]) < 100 else 'L'
            lab = f'{_LABEL.get(pm.names[jph], pm.names[jph])} + {fl}'
        else:
            lab = f'{_LABEL.get(pm.names[i])} + {_LABEL.get(pm.names[jph])}'
        w = np.array([wmap.get(r, 1.0) for r in pts[:, 0]])
        out.append((lab, ix(_wmedian(fx(pts[:, 0]), w)), _wmedian(pts[:, 1], w), 'two-phase'))
    return out


# ---- property of the stable phase ----------------------------------------------
MAP_PROPS = ('rho', 'P', 'G', 'S', 'U', 'H', 'A', 'Cp', 'Cv', 'Kt', 'Kp', 'Ks', 'alpha', 'vel',
             'Js', 'gamma_Gruneisen')
ICE_ONLY_PROPS = ('shear', 'Vp', 'Vs')                  # NaN in the fluid


def property_map(diag, *props, path=defpath):
    """Property of the stable phase over a phase diagram.

    Each grid point carries the property of the phase stable there, so the
    maps jump across the phase boundaries (density at melting and boiling...).

    * DiagramPT: every phase is evaluated on the diagram's (P, T) grid.
    * DiagramRhoT: the fluid is evaluated directly at (rho, T) (Helmholtz);
      an ice field is interpolated in density along each isotherm of the
      underlying (P, T) map, where the ice density rises with P.  Two-phase
      regions are NaN.

    Below the fluid's lowest temperature the vapour is its dilute ideal-gas
    part, as in phase_map (flagged in `ideal_gas`).

    :param diag:  DiagramPT or DiagramRhoT
    :param props: names from MAP_PROPS, plus ICE_ONLY_PROPS (NaN in the
                  fluid); none = MAP_PROPS
    :return:      PropertyMap
    """
    props = tuple(props) or MAP_PROPS
    bad = set(props) - set(MAP_PROPS) - set(ICE_ONLY_PROPS)
    if bad:
        raise ValueError('unsupported property name(s): ' + ', '.join(sorted(bad)))
    if isinstance(diag, DiagramPT):
        pm = diag.pm
        vals, ig = _stable_props_PT(pm, props, range(len(pm.names)), path)
        return PropertyMap(values=vals, phase=pm.stable.copy(), ideal_gas=ig, coords='PT',
                           names=list(pm.names))
    if isinstance(diag, DiagramRhoT):
        return _property_map_rhoT(diag, props, path)
    raise TypeError('diag must be a DiagramPT or a DiagramRhoT')


def _stable_props_PT(pm, props, phases, path):
    """{prop: (nP, nT)} of the stable phase (only `phases`), and the ideal-gas mask."""
    shape = pm.stable.shape
    vals = {p: np.full(shape, np.nan) for p in props}
    ig = np.zeros(shape, bool)
    Pm, Tm = np.meshgrid(pm.P, pm.T, indexing='ij')
    if 'P' in vals:
        vals['P'] = np.where(np.isin(pm.stable, list(phases)), Pm, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for k in phases:
            sel = pm.stable == k
            if not sel.any():
                continue
            name = pm.names[k]
            want = [p for p in props if p != 'P' and (k > 0 or p not in ICE_ONLY_PROPS)]
            if want:
                o = getProp(_grid(pm.P, pm.T), name, path, *want)
                for p in want:
                    a = np.broadcast_to(np.asarray(getattr(o, p, np.nan), float), shape)
                    vals[p][sel] = a[sel]
            if k == 0:
                sp = _load_spline(path, name)
                if eh.is_psi(sp):
                    lo = sel & (Tm < eh.domain(sp)[1][0]) & (Pm < 1e-4)   # as phase_map
                    if lo.any():
                        ig |= lo
                        want = [p for p in props if p in eh.IDEAL_GAS_PROPS and p != 'P']
                        if want:
                            o = eh.ideal_gas_props(sp, Pm[lo], Tm[lo], *want)
                            for p in want:
                                vals[p][lo] = getattr(o, p)
    return vals, ig


def _property_map_rhoT(d, props, path):
    pm, fluid = d.pm, d.pm.names[0]
    shape = d.field.shape
    vals = {p: np.full(shape, np.nan) for p in props}
    ig = np.zeros(shape, bool)
    Rm, Tm = np.meshgrid(d.rho, d.T, indexing='ij')
    # the fluid: directly at (rho, T)
    sel = d.field == 0
    if sel.any():
        if 'rho' in vals:
            vals['rho'][sel] = Rm[sel]
        want = [p for p in props if p not in ICE_ONLY_PROPS and p != 'rho']
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            if want:
                o = getProp(_grid(d.rho, d.T), fluid, path, *want, rhoT=True)
                for p in want:
                    a = np.broadcast_to(np.asarray(getattr(o, p, np.nan), float), shape)
                    vals[p][sel] = a[sel]
            sp = _load_spline(path, fluid)
            if eh.is_psi(sp):
                lo = sel & (Tm < eh.domain(sp)[1][0])          # dilute-vapour extension
                if lo.any():
                    ig |= lo
                    want = [p for p in props if p in eh.IDEAL_GAS_PROPS and p != 'rho']
                    if want:
                        o = eh.ideal_gas_props(sp, Rm[lo], Tm[lo], *want, rhoT=True)
                        for p in want:
                            vals[p][lo] = getattr(o, p)
    # the ices: interpolate in density along each isotherm of the (P, T) map
    ices = [k for k in range(1, len(pm.names)) if (d.field == k).any()]
    if ices:
        need = tuple(p for p in props if p not in ('rho', 'P'))
        src, _ = _stable_props_PT(pm, need, ices, path)
        for k in ices:
            for j in np.flatnonzero((d.field == k).any(axis=0)):
                col = pm.stable[:, j] == k
                tgt = d.field[:, j] == k
                if not col.any():
                    continue
                r = np.maximum.accumulate(pm.rho[k][col, j])
                x = d.rho[tgt]
                if 'rho' in vals:
                    vals['rho'][tgt, j] = x
                if 'P' in vals:
                    vals['P'][tgt, j] = np.exp(np.interp(x, r, np.log(pm.P[col])))
                for p in need:
                    vals[p][tgt, j] = np.interp(x, r, src[p][col, j])
    return PropertyMap(values=vals, phase=d.field.copy(), ideal_gas=ig, coords='rhoT',
                       names=list(pm.names))


# ---- the diagrams as matplotlib figures ------------------------------------------
def wpd_PT(ax=None, P=(1e-8, 1e5), T=(150.0, 1800.0), nP=500, nT=420, fluid='water3',
           ices=ICES, path=defpath, return_map=False):
    """Full H2O phase diagram in (P, T): vapour, liquid, supercritical fluid, ices.

    Draws phase_diagram_PT: stability fields by Gibbs-energy minimisation on a
    log-P grid, the boundaries between neighbouring stable phases, the
    saturation curve and critical point, and the triple points.

    :return: matplotlib Figure (and the PhaseMap if return_map=True)
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    d = phase_diagram_PT(P, T, nP, nT, fluid, ices, path)
    pm = d.pm
    if ax is None:
        fig, ax = plt.subplots(figsize=(11, 7.5))
    else:
        fig = ax.figure
    cmap = ListedColormap([_tint(c) for c in _COLORS[:len(pm.names)]])
    S = np.ma.masked_less(pm.stable, 0).astype(float)
    ax.pcolormesh(pm.P, pm.T, S.T, cmap=cmap, vmin=-0.5, vmax=len(pm.names) - 0.5, shading='auto',
                  rasterized=True)
    _hatch_missing(ax, pm.P, pm.T, d.missing)
    for i, j, Pl, Tl in d.boundaries:
        ax.plot(Pl, Tl, '-', color='#0b0b0b', lw=1.1)
    # vapour-liquid: saturation curve and critical point
    if d.saturation is not None:
        ax.plot(d.saturation.P, d.saturation.T, '-', color='#0b0b0b', lw=1.1)
    ax.plot(_PC, _TC, 'o', color='#0b0b0b', ms=6, mfc='white', mew=1.5, zorder=5)
    ax.annotate('critical point', (_PC, _TC), xytext=(8, -12), textcoords='offset points', fontsize=9)
    for tp in d.triple_points:
        ax.plot(tp['P'], tp['T'], 'o', color='#0b0b0b', ms=3.5, zorder=6)
    ax.set_xscale('log')
    ax.set_xlim(pm.P[0], pm.P[-1]); ax.set_ylim(pm.T[0], pm.T[-1])
    ax.set_xlabel('Pressure (MPa)'); ax.set_ylabel('Temperature (K)')
    ax.set_title(f'H$_2$O phase diagram — fluid: {fluid}, ices: ' + ', '.join(_LABEL[i] for i in ices),
                 loc='left', fontsize=11)
    for text, Px, Tx, kind in d.labels:
        if kind == 'fluid':
            ax.text(Px, Tx, text, ha='center', va='center', fontsize=11, fontweight='bold',
                    color='#123c70')
        else:
            ax.text(Px, Tx, text, ha='center', va='center', fontsize=10, fontweight='bold',
                    color='#0b0b0b')
    _style(ax)
    return (fig, pm) if return_map else fig


def _hatch_missing(ax, X, Y, missing):
    """Hatch and label the grid cells where no phase is available."""
    if not missing.any():
        return
    ax.contourf(X, Y, missing.T.astype(float), levels=[0.5, 1.5], colors='none', hatches=['////'])
    ax.contour(X, Y, missing.T.astype(float), levels=[0.5], colors='#8c8b87', linewidths=0.8)
    lX, YY = np.meshgrid(np.log10(X), Y, indexing='ij')
    ax.text(10 ** np.median(lX[missing]), np.median(YY[missing]), 'not modelled\n(ice VII/X field)',
            ha='center', va='center', fontsize=9, color='#52514e',
            bbox=dict(fc='white', ec='none', alpha=0.85, pad=2))


def _sat_T(Tlo, n=80):
    """Temperatures for the saturation curve: uniform, plus clustering toward Tc."""
    return np.unique(np.r_[np.linspace(Tlo, _TC - 1.0, n), _TC - np.geomspace(0.05, 30.0, 40)])


def pm_stable_fluid(pm, P, T):
    """True where the fluid (index 0) is the stable phase at the (P, T) points."""
    iP = np.clip(np.searchsorted(pm.P, P), 0, pm.P.size - 1)
    iT = np.clip(np.searchsorted(pm.T, T), 0, pm.T.size - 1)
    return pm.stable[iP, iT] == 0


def _fluid_kind(pm):
    """0 vapour, 1 liquid, 2 supercritical for fluid points; -1 elsewhere."""
    Tm = np.broadcast_to(pm.T[None, :], pm.stable.shape)
    k = np.where(pm.rho[0] >= _RHOC, 1, 0)
    k = np.where(Tm >= _TC, 2, k)
    return np.where(pm.stable == 0, k, -1)


def wpd_rhoT(ax=None, rho=(1e-7, 4e3), T=(150.0, 1800.0), P=(1e-10, 1e5), nP=1200, nT=420,
             nrho=800, xscale='log', fluid='water3', ices=ICES, path=defpath, return_map=False):
    """Full H2O phase diagram in (rho, T).

    Draws phase_diagram_rhoT: the stability fields re-drawn in density, the
    two-phase regions (grey) bounded by the coexisting densities, the
    saturation dome and critical point, and the three-phase tie lines.

    :param xscale: 'log' (default, shows the vapour) or 'linear'
    :return: matplotlib Figure (and the PhaseMap if return_map=True)
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    d = phase_diagram_rhoT(rho, T, P, nP, nT, nrho, xscale, fluid, ices, path)
    pm, rq, Tg, n = d.pm, d.rho, d.T, d.two_phase
    if ax is None:
        fig, ax = plt.subplots(figsize=(11, 7.5))
    else:
        fig = ax.figure
    cmap = ListedColormap([_tint(c) for c in _COLORS[:n]] + [(0.86, 0.86, 0.84)])
    ax.pcolormesh(rq, Tg, np.ma.masked_less(d.field, 0).T, cmap=cmap, vmin=-0.5, vmax=n + 0.5,
                  shading='auto', rasterized=True)
    for i, j, ri, rj, Tl in d.coexistence:
        ax.plot(ri, Tl, '-', color='#0b0b0b', lw=0.8)
        ax.plot(rj, Tl, '-', color='#0b0b0b', lw=0.8)
    if d.saturation is not None:
        ax.plot(d.saturation.rho_A, d.saturation.T, '-', color='#0b0b0b', lw=1.1)
        ax.plot(d.saturation.rho_B, d.saturation.T, '-', color='#0b0b0b', lw=1.1)
    ax.plot(_RHOC, _TC, 'o', color='#0b0b0b', ms=6, mfc='white', mew=1.5, zorder=5)
    ax.annotate('critical point', (_RHOC, _TC), xytext=(-14, 10), textcoords='offset points',
                fontsize=9, ha='right')
    # three-phase (triple-point) tie lines: the separations between two-phase regions
    for tp in d.tie_lines:
        r = np.array(tp['rho'])
        ax.plot([r.min(), r.max()], [tp['T'], tp['T']], '-', color='#0b0b0b', lw=1.0, zorder=4)
        ax.plot(r, np.full(3, tp['T']), '|', color='#0b0b0b', ms=7, mew=1.2, zorder=4)
    ax.set_xscale(xscale)
    ax.set_xlim(rq[0], rq[-1]); ax.set_ylim(Tg[0], Tg[-1])
    ax.set_xlabel('Density (kg/m³)'); ax.set_ylabel('Temperature (K)')
    ax.set_title(f'H$_2$O phase diagram in density — fluid: {fluid}; grey: two-phase regions',
                 loc='left', fontsize=11)
    for text, x, Tx, kind in d.labels:
        if kind == 'fluid':
            ax.text(x, Tx, text, ha='center', va='center', fontsize=11, fontweight='bold',
                    color='#123c70')
        elif kind == 'ice':
            ax.text(x, Tx, text, ha='center', va='center', fontsize=10, fontweight='bold')
        else:
            ax.text(x, Tx, text, ha='center', va='center', fontsize=8.5, color='#52514e',
                    style='italic')
    _style(ax)
    return (fig, pm) if return_map else fig


def _rho_on(pm, ph, P, T, fluid, path):
    """Density of phase `ph` at scattered (P, T) points on a boundary."""
    pts = np.empty(P.size, dtype=object)
    for q in range(P.size):
        pts[q] = (float(P[q]), float(T[q]))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        if ph == 0:
            r = np.asarray(getProp(pts, fluid, path, 'rho').rho, float)
            lo = ~np.isfinite(r)
            if lo.any():                                  # dilute-vapour extension
                sp = _load_spline(path, fluid)
                r[lo] = eh.ideal_gas(sp, P[lo], T[lo])[1]
            return r
        return np.asarray(getProp(pts, pm.names[ph], path, 'rho').rho, float)
