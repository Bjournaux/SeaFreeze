from collections import namedtuple
from itertools import repeat
import types
import warnings
import logging
import os.path as op
import numpy as np
from mlbspline import load
from lbftd import evalGibbs as eg
from lbftd.statevars import iP, iT, iM

log = logging.getLogger('seafreeze')
stream = logging.StreamHandler()
stream.setFormatter(logging.Formatter('[SeaFreeze %(levelname)s] %(message)s'))
defpath = op.join(op.dirname(op.abspath(__file__)), 'splines')

# Map from material code → (subfolder, filename)
# Mirrors the Matlab sf_load_spline map; all individual files store the spline as 'sp'.
_SPLINE_MAP = {
    'Ih':              ('ice_Ih',              'ice_Ih.mat'),
    'II':              ('ice_II',              'ice_II.mat'),
    'III':             ('ice_III',             'ice_III.mat'),
    'V':               ('ice_V',               'ice_V.mat'),
    'VI':              ('ice_VI',              'ice_VI.mat'),
    'VII_X_French':    ('ice_VII_X_French',    'ice_VII_X_French.mat'),
    'water1':          ('water_Bollengier',     'water_Bollengier.mat'),
    'water2':          ('water_Brown',          'water_Brown.mat'),
    'water_IAPWS95':   ('water_IAPWS95',        'water_IAPWS95.mat'),
    'NaClaq':          None,  # stitched LP+HP — handled separately
    'NaClaq_LP':       ('NaCl_aq_LP_2026',     'NaCl_aq_LP_2026.mat'),
    'NaClaq_HP':       ('NaCl_aq_HP_2026',     'NaCl_aq_HP_2026.mat'),
    'NaClaq_5GPa_2024':('NaCl_aq_Brown2024',   'NaCl_aq_Brown2024.mat'),
}

# ---------------------------------------------------------------------------
# Property classification — mirrors Matlab SF_getprop property lists
# ---------------------------------------------------------------------------
# Properties we compute here (not by lbftd) from existing lbftd outputs.
# Base props: available for all phases.
_DERIVED_BASE = frozenset({'Js', 'gamma_Gruneisen', 'P', 'T'})
# NaClaq-only derived props.
_DERIVED_NACL = frozenset({'m', 'xs', 'xw', 'f', 'Vw'})
# lbftd prerequisite props required for each SF-derived prop.
_PREREQS = {
    'Js':              frozenset({'alpha', 'rho', 'Cp'}),
    'gamma_Gruneisen': frozenset({'alpha', 'Kt',  'rho', 'Cv'}),
    'Vw':              frozenset({'muw'}),
}

# ---------------------------------------------------------------------------
# Constants (module-level; used in formulas below)
# ---------------------------------------------------------------------------
mH2O_kgmol  = 18.01528e-3   # kg/mol  — water
mNaCl_kgmol = 58.44e-3      # kg/mol  — NaCl
# NaClaq LP+HP stitching boundaries (MPa) — mirrors SF_NaCl_stitch.m
_STITCH_P_LO = 499.9    # LP upper knot boundary / start of cosine taper
_STITCH_P_HI = 1001.0   # HP lower knot boundary / end of cosine taper
# P-independent fields: depend only on molality, copied from either spline
_P_INDEP_FIELDS = frozenset({'xs', 'xw'})
# Coordinate-echo fields: set by _compute_derived, not assembled
_COORD_FIELDS = frozenset({'P', 'T', 'm'})
# Complete set of lbftd props to request in all-props mode.
# V, gam, Gex are lbftd-only and not part of the Matlab output — never requested.
_LBFTD_ALL_WATER = ('G', 'S', 'U', 'H', 'A', 'rho', 'Cp', 'Cv',
                    'Kt', 'Kp', 'Ks', 'alpha', 'vel')
_LBFTD_ALL_NACL  = _LBFTD_ALL_WATER + (
                    'mus', 'muw', 'Vm', 'Cpm', 'phi', 'aw',
                    'Va', 'Vex', 'Cpa')

# mol of water per kg of water = 1 kg / (18.01528e-3 kg/mol) ≈ 55.508 mol/kg
_nw = 1.0 / mH2O_kgmol      # mol water / kg water ≈ 55.508

PhaseDesc = namedtuple('PhaseDesc', 'shear_mod_parms phase_num MW nu cutoff')
phases = {
    # Ice phases — Journaux et al. 2020 / Feistel & Wagner 2006
    "Ih":           PhaseDesc([3.1, -0.00462, 0, -0.00657, 1000, 273.15], 1, mH2O_kgmol, None, None),
    "II":           PhaseDesc([4.1,  0.0175,  0, -0.014,   1100, 273],    2, mH2O_kgmol, None, None),
    "III":          PhaseDesc([2.57, 0.0175,  0, -0.014,   1100, 273],    3, mH2O_kgmol, None, None),
    "V":            PhaseDesc([2.57, 0.0175,  0, -0.014,   1100, 273],    5, mH2O_kgmol, None, None),
    "VI":           PhaseDesc([2.57, 0.0175,  0, -0.014,   1100, 273],    6, mH2O_kgmol, None, None),
    "VII_X_French": PhaseDesc([10,   0.0033,  0.000048, -0.014, 1300, 273], 7, mH2O_kgmol, None, None),
    # Pure water — liquid phases
    "water1":        PhaseDesc(None, 0,      mH2O_kgmol, None, None),  # Bollengier et al. 2019 (≤500 K, ≤2300 MPa)
    "water2":        PhaseDesc(None, np.nan, mH2O_kgmol, None, None),  # Brown 2018 (up to 100 GPa)
    "water_IAPWS95": PhaseDesc(None, np.nan, mH2O_kgmol, None, None),  # IAPWS95; Wagner & Pruss 2002
    # Aqueous NaCl
    "NaClaq":          PhaseDesc(None, 0, mNaCl_kgmol, 2, 0.0002),  # stitched LP+HP 2026 (recommended)
    "NaClaq_LP":       PhaseDesc(None, 0, mNaCl_kgmol, 2, 0.0002),  # 2026 low-P  spline only
    "NaClaq_HP":       PhaseDesc(None, 0, mNaCl_kgmol, 2, 0.0002),  # 2026 high-P spline only
    "NaClaq_5GPa_2024":PhaseDesc(None, 0, mNaCl_kgmol, 2, 0.0002),  # Brown 2024 legacy spline
}
max_phase_num = int(np.nanmax([p.phase_num for p in phases.values()]))

# Build phase_num → material code map; exclude NaN phase_nums and keep only
# the first entry for phase_num 0 (water1 is canonical; NaClaq variants share 0).
phasenum2phaseDict = {}
for k, v in phases.items():
    if not np.isnan(v.phase_num) and v.phase_num not in phasenum2phaseDict:
        phasenum2phaseDict[int(v.phase_num)] = k


# ---------------------------------------------------------------------------
# Spline loading
# ---------------------------------------------------------------------------
def _load_spline(splines_dir, material):
    """Load the Gibbs LBF spline for *material* from the per-spline folder structure.

    Mirrors the Matlab ``sf_load_spline`` function.  All individual spline files
    store the spline struct under the variable name ``sp``.

    :param splines_dir: path to the ``splines/`` directory (= defpath)
    :param material:    material code (key of _SPLINE_MAP)
    :return:            spline dict as returned by mlbspline.load.loadSpline
    """
    entry = _SPLINE_MAP.get(material)
    if entry is None:
        if material == 'NaClaq':
            raise ValueError(
                "'NaClaq' uses LP+HP stitching and cannot be loaded as a single "
                "spline.  Use getProp(PTm, 'NaClaq') which handles stitching "
                "automatically, or load 'NaClaq_LP' / 'NaClaq_HP' individually.")
        else:
            raise ValueError(f"Unknown material '{material}'. Supported: "
                             + ', '.join(k for k, v in _SPLINE_MAP.items() if v is not None))
    subfolder, filename = entry
    sp = load.loadSpline(op.join(splines_dir, subfolder, filename), 'sp')
    # Ensure NaClaq splines have metadata lbftd expects
    pd = phases.get(material)
    if pd is not None and pd.cutoff is not None:
        if 'cutoff' not in sp or sp['cutoff'] is None:
            sp['cutoff'] = pd.cutoff
        # lbftd expects MW as [MW_solvent, MW_solute]; some splines store only MW_solute
        mw = sp.get('MW')
        if mw is not None and np.isscalar(mw):
            sp['MW'] = np.array([mH2O_kgmol, float(mw)])
    return sp


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def seafreeze(PTm, phase, path=defpath, *tdvSpec):
    # this deprecation warning is being added on 2024.06.21.
    # Per Python guidelines, the function should not be removed for 2 years.
    warnings.warn('use getProp instead', DeprecationWarning, stacklevel=2)
    return getProp(PTm, phase, path, *tdvSpec)


def getProp(PTm, phase, path=defpath, *tdvSpec, verbose=False):
    """Calculates thermodynamic quantities for H2O water or ice polymorphs
    Ih, II, III, V, VI, VII/X and aqueous NaCl.

    Output properties match the Matlab ``SF_getprop`` reference:

    All phases
    ----------
    G (J/kg), S (J/K/kg), U (J/kg), H (J/kg), A (J/kg),
    rho (kg/m³), Cp (J/kg/K), Cv (J/kg/K),
    Kt (MPa), Kp (−), Ks (MPa), alpha (1/K), vel (m/s),
    Js (K/MPa), gamma_Gruneisen (−), P (MPa echo), T (K echo)

    Solid phases additionally
    -------------------------
    shear (MPa), Vp (m/s), Vs (m/s)

    NaClaq additionally
    -------------------
    mus, muw, Va, Cpa, Vm, Cpm, phi, Vex, aw,
    m (mol/kg echo), xs, xw, f, Vw (cm³/mol)

    NOTE:  The authors recommend 'water1' for 200–355 K up to 2300 MPa.
    The ice Gibbs parametrizations are optimized for phase-equilibrium
    calculations against 'water1'.  'water2' and 'water_IAPWS95' are
    provided for HP extension and comparison only.

    :param PTm:     P (MPa) / T (K) conditions, optional m (mol/kg) for NaClaq.
                    Scatter: 1-D numpy array of (P,T) or (P,T,m) tuples.
                    Grid: numpy array([P_vec, T_vec]) or ([P, T, m_vec]).
    :param phase:   Material code — key of the ``phases`` dict.
    :param path:    Path to the ``splines/`` directory (default: package splines).
    :param tdvSpec: Optional property names to compute; default = all supported.
    :return:        Object with computed properties as named attributes.
    """
    lbftd_log = logging.getLogger('lbftd')
    if not verbose:
        lbftd_log.setLevel(logging.CRITICAL)
    try:
        try:
            phasedesc = phases[phase]
        except KeyError:
            raise ValueError('The specified phase is not recognized.  Supported phases are ' +
                             ', '.join(phases.keys()) + '.')

        isscatter = _is_scatter(PTm)
        want_set  = set(tdvSpec)
        want_all  = len(want_set) == 0
        is_nacl   = phase.startswith('NaClaq')
        is_stitched = (phase == 'NaClaq')

        # ---- Determine which derived props are relevant for this phase ----------
        derived_known = _DERIVED_BASE | (_DERIVED_NACL if is_nacl else set())

        # ---- Build the set of props to request from lbftd ----------------------
        # V, gam, Gex are lbftd-internal and never exposed (matching Matlab).
        if want_all:
            lbftd_wants = _LBFTD_ALL_NACL if is_nacl else _LBFTD_ALL_WATER
        else:
            lbftd_wants = tuple(want_set - derived_known - {'shear', 'Vp', 'Vs'})
            extra = set()
            if 'Js'              in want_set: extra |= _PREREQS['Js']
            if 'gamma_Gruneisen' in want_set: extra |= _PREREQS['gamma_Gruneisen']
            if 'Vw'              in want_set: extra |= _PREREQS['Vw']
            if want_set & {'shear', 'Vp', 'Vs'} and phasedesc.shear_mod_parms:
                extra |= {'rho', 'Ks'}
            lbftd_wants = tuple(set(lbftd_wants) | extra)

        # ---- Evaluate — either stitched LP+HP or single spline -----------------
        if is_stitched:
            sp = None  # no single spline for stitched mode
            props = _nacl_stitch(PTm, isscatter, path, *lbftd_wants)
        else:
            sp = _load_spline(path, phase)
            raw = _get_tdvs(sp, PTm, isscatter, *lbftd_wants)
            # Use instance __dict__ to avoid lbftd class-level attribute fallback.
            props = dict(vars(raw))

        # ---- Shear / Vp / Vs (solid phases only) --------------------------------
        if phasedesc.shear_mod_parms and (want_all or want_set & {'shear', 'Vp', 'Vs'}):
            Ks  = np.asarray(props['Ks'])
            rho = np.asarray(props['rho'])
            smg = _get_shear_mod_GPa(phasedesc.shear_mod_parms, rho, _get_T(PTm, isscatter))
            if want_all or 'shear' in want_set: props['shear'] = 1e3 * smg
            if want_all or 'Vp'    in want_set: props['Vp']    = _get_Vp(smg, rho, Ks)
            if want_all or 'Vs'    in want_set: props['Vs']    = _get_Vs(smg, rho)
            if not want_all:
                if 'Ks'  not in want_set: props.pop('Ks', None)
                if 'rho' not in want_set: props.pop('rho', None)

        # ---- Matlab-parity derived properties -----------------------------------
        _compute_derived(props, PTm, isscatter, sp, phase, phasedesc, want_set, path)

        # ---- Strip prerequisite props that were only added internally -----------
        if not want_all:
            for p in ('alpha', 'Cp', 'Cv', 'Kt', 'muw'):
                if p not in want_set:
                    props.pop(p, None)

        return types.SimpleNamespace(**props)
    finally:
        if not verbose:
            lbftd_log.setLevel(logging.WARNING)


def whichphase(PTm, solute='water1', path=defpath):
    """Determines the most likely phase of water at each pressure/temperature.

    :param PTm:     P (MPa) / T (K) conditions (and optional m for NaClaq).
                    Scatter: 1-D numpy array of (P,T) or (P,T,m) tuples.
                    Grid: numpy array([P_vec, T_vec]) or ([P, T, m_vec]).
    :param solute:  An optional dissolved solute in the liquid phase.
                    The default is pure water (water1).
    :param path:    Path to the ``splines/`` directory.
    :return:        numpy.ndarray with the stable phase index at each point.
    """
    is_stitched_solute = (solute == 'NaClaq')
    isscatter = _is_scatter(PTm)
    # NaClaq uses stitched LP+HP evaluation and cannot be loaded as a single spline.
    # Load only solid-ice splines (phase_num > 0); liquid is handled below.
    if is_stitched_solute:
        phase_sp = {v.phase_num: _load_spline(path, pcomp) for pcomp, v in phases.items() if
                    v.phase_num > 0}
    else:
        phase_sp = {v.phase_num: _load_spline(path, pcomp) for pcomp, v in phases.items() if
                    v.phase_num > 0 or pcomp == solute}
    ptsh = (PTm.size,) if isscatter else (PTm[iP].size, PTm[iT].size)
    comp = np.full(ptsh + (max_phase_num + 1,), np.nan)
    for p in phase_sp.keys():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sl = tuple(repeat(slice(None), 1 if isscatter else 2)) + (p,)
            sp = phase_sp[p]
            if p == 0:
                if 'water1' in solute:
                    phase_sp[p]['MW'] = phases[solute].MW
                    tdvs = _get_tdvs(sp, _get_PT(PTm, isscatter), isscatter, 'G').G * phase_sp[p]['MW']
                else:
                    phase_sp[p]['nu'] = phases[phasenum2phase(p)].nu
                    tdvs = _get_tdvs(sp, PTm, isscatter, 'G', 'muw').muw
            else:
                phase_sp[p]['MW'] = phases[phasenum2phase(p)].MW
                tdvs = _get_tdvs(sp, _get_PT(PTm, isscatter), isscatter, 'G').G * phase_sp[p]['MW']
            comp[sl] = np.squeeze(tdvs)
    # Stitched NaClaq liquid: evaluate muw via getProp (handles LP+HP blending)
    if is_stitched_solute:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            sl = tuple(repeat(slice(None), 1 if isscatter else 2)) + (0,)
            comp[sl] = np.squeeze(getProp(PTm, 'NaClaq', path, 'muw').muw)
    all_nan_sl = np.all(np.isnan(comp), -1)
    out = np.full(ptsh, np.nan)
    out[~all_nan_sl] = np.nanargmin(comp[~all_nan_sl], -1)
    return out


def phasenum2phase(phaseInt, liqComp='water1'):
    """Convert an integer phase index to a material string.

    :param phaseInt: Integer representing the desired ice phase or liquid.
    :param liqComp: String to return for the liquid phase (phase_num = 0).
    :return: Material string compatible with getProp/whichphase.
    """
    if phaseInt == 0:
        return liqComp
    elif np.isnan(phaseInt):
        log.warning('NaN input to phasenum2phase. Defaulting to liquid phase.')
        return liqComp
    else:
        return phasenum2phaseDict[phaseInt]


# ---------------------------------------------------------------------------
# NaClaq LP+HP stitching (mirrors Matlab SF_NaCl_stitch.m)
# ---------------------------------------------------------------------------
def _nacl_stitch(PTm, isscatter, path, *tdvSpec):
    """Evaluate NaClaq properties by stitching LP and HP splines.

    In the overlap zone [499.9, 1001] MPa each property is blended with a
    C-infinity cosine taper:  w(P) = 0.5*(1 + cos(pi*(P-P_LO)/(P_HI-P_LO)))
    Outside: LP only (P <= P_LO) or HP only (P >= P_HI).

    Returns a plain dict of property arrays (same format as getProp's internal
    ``props`` dict).
    """
    sp_lp = _load_spline(path, 'NaClaq_LP')
    sp_hp = _load_spline(path, 'NaClaq_HP')
    # Ensure cutoff is set (lbftd needs it for apparent-property evaluation)
    cutoff = phases['NaClaq'].cutoff
    if 'cutoff' not in sp_lp or sp_lp['cutoff'] is None:
        sp_lp['cutoff'] = cutoff
    if 'cutoff' not in sp_hp or sp_hp['cutoff'] is None:
        sp_hp['cutoff'] = cutoff

    P_arr = np.asarray(_get_P(PTm, isscatter), dtype=float).ravel()
    nP = P_arr.size

    # --- Cosine taper weights ---
    w_P = np.zeros(nP)
    w_P[P_arr <= _STITCH_P_LO] = 1.0
    # P >= P_HI already 0
    in_taper = (P_arr > _STITCH_P_LO) & (P_arr < _STITCH_P_HI)
    w_P[in_taper] = 0.5 * (1.0 + np.cos(
        np.pi * (P_arr[in_taper] - _STITCH_P_LO) / (_STITCH_P_HI - _STITCH_P_LO)))

    # --- Logical masks: which P values each spline is evaluated at ---
    need_lp = P_arr <= _STITCH_P_HI
    need_hp = P_arr >= _STITCH_P_LO

    # --- Build sub-inputs ---
    if isscatter:
        sub_lp = PTm[need_lp] if np.any(need_lp) else None
        sub_hp = PTm[need_hp] if np.any(need_hp) else None
    else:
        P_full = np.asarray(PTm[0], dtype=float).ravel()
        T_vec  = PTm[1]
        m_vec  = PTm[2]
        sub_lp = np.array([P_full[need_lp], T_vec, m_vec], dtype=object) if np.any(need_lp) else None
        sub_hp = np.array([P_full[need_hp], T_vec, m_vec], dtype=object) if np.any(need_hp) else None

    # --- Evaluate each sub-spline ---
    # allowExtrapolations=True mirrors Matlab fnGval behaviour: the union T
    # domain reported by phase_range includes a thin sliver (≤0.9 K) outside
    # the LP spline's lower T knot, and fnGval extrapolates there gracefully.
    fn = eg.evalSolutionGibbsScatter if isscatter else eg.evalSolutionGibbsGrid
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res_lp = dict(vars(fn(sp_lp, sub_lp, *tdvSpec, allowExtrapolations=True))) if sub_lp is not None else {}
        res_hp = dict(vars(fn(sp_hp, sub_hp, *tdvSpec, allowExtrapolations=True))) if sub_hp is not None else {}

    # --- Index mappings (P dimension only) ---
    idx_lp = np.where(need_lp)[0]
    idx_hp = np.where(need_hp)[0]

    blend_global = np.where(need_lp & need_hp)[0]
    lp_only      = np.where(need_lp & ~need_hp)[0]
    hp_only      = np.where(~need_lp & need_hp)[0]

    # Positions within LP / HP result arrays
    blend_in_lp = np.searchsorted(idx_lp, blend_global)
    blend_in_hp = np.searchsorted(idx_hp, blend_global)
    lp_only_in  = np.searchsorted(idx_lp, lp_only)
    hp_only_in  = np.searchsorted(idx_hp, hp_only)

    # Taper weights for the blend zone
    w_blend = w_P[blend_global]

    # --- Assemble output dict field by field ---
    all_fields = set(res_lp.keys()) | set(res_hp.keys())
    props = {}

    for fn in all_fields:
        if fn in _COORD_FIELDS:
            continue  # handled later by _compute_derived
        if fn in _P_INDEP_FIELDS:
            # xs, xw depend only on molality — copy from either spline
            props[fn] = res_lp.get(fn, res_hp.get(fn))
            continue

        arr_lp = res_lp.get(fn)
        arr_hp = res_hp.get(fn)

        if isscatter:
            props[fn] = _assemble_scatter(
                arr_lp, arr_hp, nP,
                lp_only, lp_only_in, hp_only, hp_only_in,
                blend_global, blend_in_lp, blend_in_hp, w_blend)
        else:
            props[fn] = _assemble_grid(
                arr_lp, arr_hp, nP,
                lp_only, lp_only_in, hp_only, hp_only_in,
                blend_global, blend_in_lp, blend_in_hp, w_blend)

    return props


def _assemble_grid(arr_lp, arr_hp, nP,
                   lp_only, lp_only_in, hp_only, hp_only_in,
                   blend_idx, blend_in_lp, blend_in_hp, w_blend):
    """Assemble an (nP, nT, ...) array from LP and HP sub-arrays with blending.

    The first axis is always the P axis.  Trailing dimensions (T, m) are
    preserved.  Mirrors Matlab ``assemble_grid`` in ``SF_NaCl_stitch.m``.
    """
    # Infer trailing shape from whichever result is available
    ref = arr_lp if arr_lp is not None else arr_hp
    ref = np.asarray(ref)
    tail = ref.shape[1:]
    ntail = int(np.prod(tail)) if tail else 1

    out2 = np.zeros((nP, ntail))

    if lp_only.size and arr_lp is not None:
        lp2 = np.asarray(arr_lp).reshape(arr_lp.shape[0], ntail) if ntail > 1 else np.asarray(arr_lp).reshape(-1, 1)
        out2[lp_only, :] = lp2[lp_only_in, :]
    if hp_only.size and arr_hp is not None:
        hp2 = np.asarray(arr_hp).reshape(arr_hp.shape[0], ntail) if ntail > 1 else np.asarray(arr_hp).reshape(-1, 1)
        out2[hp_only, :] = hp2[hp_only_in, :]
    if blend_idx.size and arr_lp is not None and arr_hp is not None:
        lp2 = np.asarray(arr_lp).reshape(arr_lp.shape[0], ntail) if ntail > 1 else np.asarray(arr_lp).reshape(-1, 1)
        hp2 = np.asarray(arr_hp).reshape(arr_hp.shape[0], ntail) if ntail > 1 else np.asarray(arr_hp).reshape(-1, 1)
        w = w_blend[:, np.newaxis]
        out2[blend_idx, :] = w * lp2[blend_in_lp, :] + (1.0 - w) * hp2[blend_in_hp, :]

    return out2.reshape((nP,) + tail)


def _assemble_scatter(arr_lp, arr_hp, nP,
                      lp_only, lp_only_in, hp_only, hp_only_in,
                      blend_idx, blend_in_lp, blend_in_hp, w_blend):
    """Assemble an (nP,) vector from LP and HP sub-vectors with blending.

    Mirrors Matlab ``assemble_scatter`` in ``SF_NaCl_stitch.m``.
    """
    out = np.zeros(nP)

    if lp_only.size and arr_lp is not None:
        a_lp = np.asarray(arr_lp).ravel()
        out[lp_only] = a_lp[lp_only_in]
    if hp_only.size and arr_hp is not None:
        a_hp = np.asarray(arr_hp).ravel()
        out[hp_only] = a_hp[hp_only_in]
    if blend_idx.size and arr_lp is not None and arr_hp is not None:
        a_lp = np.asarray(arr_lp).ravel()
        a_hp = np.asarray(arr_hp).ravel()
        out[blend_idx] = w_blend * a_lp[blend_in_lp] + (1.0 - w_blend) * a_hp[blend_in_hp]

    return out


# ---------------------------------------------------------------------------
# Derived-property computation (Matlab parity)
# ---------------------------------------------------------------------------
def _compute_derived(props, PTm, isscatter, sp, phase, phasedesc, want_set, path=None):
    """Attach Matlab-parity derived properties to the *props* dict in-place.

    :param props:    Plain dict built from lbftd output (modified in-place).
    :param sp:       Spline dict (None for stitched 'NaClaq').
    :param path:     Splines directory (needed for stitched Vw).
    :param want_set: set of requested property names; empty = compute all.
    """
    want_all = len(want_set) == 0

    def want(name):
        return want_all or name in want_set

    T_arr = _get_T(PTm, isscatter)
    P_arr = _get_P(PTm, isscatter)

    # For grid mode, reshape T and P so they broadcast correctly against the
    # nD prop arrays returned by lbftd.  For a 2D grid (nP×nT), T has shape
    # (nT,) which already broadcasts from the right.  For a 3D NaClaq grid
    # (nP×nT×nm), (nT,) would align with the nm axis — reshape to (1,nT,1).
    if not isscatter and props:
        sample_ndim = np.asarray(next(iter(props.values()))).ndim
        if sample_ndim == 3:          # NaClaq 3-D grid
            T_arr = np.asarray(T_arr).reshape(1, -1, 1)
            P_arr = np.asarray(P_arr).reshape(-1, 1, 1)
        # 2D case: T (nT,) already broadcasts against (nP, nT) — no change needed

    # --- Present for all phases ---
    if want('P'): props['P'] = P_arr
    if want('T'): props['T'] = T_arr

    # Js = dT/dP|adiabatic = T·α / (ρ·Cp) × 1e6   [K/MPa]
    # Derivation: Js = −(∂²G/∂P∂T) / (∂²G/∂T²) = α·d1P / (Cp/T)
    #   with d1P = 1e6/ρ → Js = α·T·1e6 / (ρ·Cp)
    if want('Js'):
        T_for_Js = props['T'] if 'T' in props else T_arr
        props['Js'] = (np.asarray(T_for_Js) *
                       np.asarray(props['alpha']) /
                       (np.asarray(props['rho']) * np.asarray(props['Cp'])) * 1e6)

    # γ_Grüneisen = α·Kt·1e6 / (ρ·Cv)   [dimensionless]
    # (1e6 converts Kt from MPa to Pa so units cancel)
    if want('gamma_Gruneisen'):
        props['gamma_Gruneisen'] = (1e6 * np.asarray(props['alpha']) *
                                    np.asarray(props['Kt']) /
                                    (np.asarray(props['rho']) * np.asarray(props['Cv'])))

    # --- NaClaq-only ---
    is_nacl = phase.startswith('NaClaq')
    if not is_nacl:
        return

    m_arr = _get_m(PTm, isscatter)
    # For 3D grid (nP×nT×nm) m_arr has shape (nm,); reshape to (1,1,nm).
    if not isscatter:
        m_arr = np.asarray(m_arr).reshape(1, 1, -1)
    nu    = phasedesc.nu        # number of ions per formula unit (NaCl → 2)
    M     = phasedesc.MW        # solute molar mass [kg/mol]

    if want('m'):  props['m']  = m_arr
    if want('xs'): props['xs'] = nu * m_arr / (nu * m_arr + _nw)
    if want('xw'): props['xw'] = 1.0 - nu * m_arr / (nu * m_arr + _nw)
    if want('f'):  props['f']  = 1.0 + M * m_arr   # kg-solution / kg-water

    # Vw: partial molar volume of water [cm³/mol]
    # = ∂muw/∂P  (central-difference, δP = 0.1 MPa)
    if want('Vw'):
        if phase == 'NaClaq':
            # Stitched mode: use _nacl_stitch at shifted P for correct blend-zone Vw
            props['Vw'] = _compute_Vw_stitched(PTm, isscatter, path)
        else:
            props['Vw'] = _compute_Vw(sp, PTm, isscatter)


def _compute_Vw(sp, PTm, isscatter, dP=0.1):
    """Partial molar volume of water [cm³/mol] via central-difference ∂muw/∂P.

    lbftd returns muw in J/mol, so ∂muw/∂P has units J/(mol·MPa) = cm³/mol
    (since 1 J/MPa = 1 cm³).  No molar-mass conversion is needed.
    """
    PTm_up = _offset_P(PTm, isscatter, +dP)
    PTm_dn = _offset_P(PTm, isscatter, -dP)
    fn = eg.evalSolutionGibbsScatter if isscatter else eg.evalSolutionGibbsGrid
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        muw_up = np.asarray(fn(sp, PTm_up, 'muw', allowExtrapolations=True).muw)
        muw_dn = np.asarray(fn(sp, PTm_dn, 'muw', allowExtrapolations=True).muw)
    # d(muw)/dP [J/(mol·MPa)] = [cm³/mol]
    return (muw_up - muw_dn) / (2.0 * dP)


def _compute_Vw_stitched(PTm, isscatter, path, dP=0.1):
    """Partial molar volume of water [cm³/mol] for stitched NaClaq.

    Mirrors Matlab SF_NaCl_stitch.m: compute Vw independently from each
    sub-spline via central-difference, then blend with the same cosine taper
    used for all other properties.  This avoids picking up the derivative of
    the blending-weight function that would appear if we finite-differenced
    the already-blended muw.
    """
    sp_lp = _load_spline(path, 'NaClaq_LP')
    sp_hp = _load_spline(path, 'NaClaq_HP')

    P_arr = np.asarray(_get_P(PTm, isscatter), dtype=float).ravel()
    nP = P_arr.size

    # Cosine taper weights at the actual P values
    w_P = np.zeros(nP)
    w_P[P_arr <= _STITCH_P_LO] = 1.0
    in_taper = (P_arr > _STITCH_P_LO) & (P_arr < _STITCH_P_HI)
    w_P[in_taper] = 0.5 * (1.0 + np.cos(
        np.pi * (P_arr[in_taper] - _STITCH_P_LO) / (_STITCH_P_HI - _STITCH_P_LO)))

    need_lp = P_arr <= _STITCH_P_HI
    need_hp = P_arr >= _STITCH_P_LO

    # Sub-inputs (same split as _nacl_stitch)
    if isscatter:
        sub_lp = PTm[need_lp] if np.any(need_lp) else None
        sub_hp = PTm[need_hp] if np.any(need_hp) else None
    else:
        P_full = np.asarray(PTm[0], dtype=float).ravel()
        T_vec, m_vec = PTm[1], PTm[2]
        sub_lp = np.array([P_full[need_lp], T_vec, m_vec], dtype=object) if np.any(need_lp) else None
        sub_hp = np.array([P_full[need_hp], T_vec, m_vec], dtype=object) if np.any(need_hp) else None

    # Vw from each sub-spline via central-difference
    vw_lp = _compute_Vw(sp_lp, sub_lp, isscatter, dP) if sub_lp is not None else None
    vw_hp = _compute_Vw(sp_hp, sub_hp, isscatter, dP) if sub_hp is not None else None

    # Blend with same index mappings as _nacl_stitch
    idx_lp = np.where(need_lp)[0]
    idx_hp = np.where(need_hp)[0]
    blend_global = np.where(need_lp & need_hp)[0]
    lp_only      = np.where(need_lp & ~need_hp)[0]
    hp_only      = np.where(~need_lp & need_hp)[0]
    blend_in_lp  = np.searchsorted(idx_lp, blend_global)
    blend_in_hp  = np.searchsorted(idx_hp, blend_global)
    lp_only_in   = np.searchsorted(idx_lp, lp_only)
    hp_only_in   = np.searchsorted(idx_hp, hp_only)
    w_blend      = w_P[blend_global]

    if isscatter:
        return _assemble_scatter(vw_lp, vw_hp, nP,
                                 lp_only, lp_only_in, hp_only, hp_only_in,
                                 blend_global, blend_in_lp, blend_in_hp, w_blend)
    else:
        return _assemble_grid(vw_lp, vw_hp, nP,
                              lp_only, lp_only_in, hp_only, hp_only_in,
                              blend_global, blend_in_lp, blend_in_hp, w_blend)


def _offset_P(PTm, is_scatter, delta):
    """Return a copy of PTm with all pressures shifted by *delta* MPa."""
    if is_scatter:
        out = np.empty(PTm.shape, dtype=object)
        for i, tup in enumerate(PTm):
            out[i] = (tup[0] + delta,) + tup[1:]
        return out
    else:
        shifted = PTm.copy()
        shifted[0] = PTm[0] + delta
        return shifted


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------
def _get_tdvs(sp, PTm, is_scatter, *tdvSpec):
    """Call the appropriate lbftd evalGibbs function."""
    fn = eg.evalSolutionGibbsScatter if is_scatter else eg.evalSolutionGibbsGrid
    return fn(sp, PTm, *tdvSpec, allowExtrapolations=False)



def _get_shear_mod_GPa(sm, rho, T):
    return None if sm is None else (
        sm[0] + sm[1] * (rho - sm[4]) + sm[2] * (rho - sm[4]) ** 2 + sm[3] * (T - sm[5]))


def _get_Vp(smg, rho, Ks):
    return 1e3 * np.sqrt((Ks / 1e3 + 4 / 3 * smg) / rho / 1e-3)


def _get_Vs(smg, rho):
    return 1e3 * np.sqrt(smg / rho / 1e-3)


def _is_scatter(PTm):
    return isinstance(PTm[0], tuple) or (PTm.shape == (1, 2) and np.isscalar(PTm[0]) and np.isscalar(PTm[1])) \
        or (PTm.shape == (1, 3) and np.isscalar(PTm[0]) and np.isscalar(PTm[1]) and np.isscalar(PTm[2]))


def _get_T(PTm, is_scatter):
    if is_scatter:
        if len(PTm[0]) < 3:
            return np.array([T for P, T in PTm])
        else:
            return np.array([T for P, T, m in PTm])
    else:
        return PTm[1]


def _get_P(PTm, is_scatter):
    if is_scatter:
        return np.array([tup[0] for tup in PTm])
    else:
        return PTm[0]


def _get_m(PTm, is_scatter):
    if is_scatter:
        return np.array([tup[2] for tup in PTm])
    else:
        return PTm[2]


def _get_PT(PTm, is_scatter):
    if is_scatter:
        if len(PTm[0]) < 3:
            return PTm
        else:
            out = np.empty((PTm.size,), object)
            for i in range(PTm.size):
                out[i] = (PTm[i][0], PTm[i][1])
            return out
    else:
        return PTm[:2]


def _del(obj, attr):
    """Remove *attr* from *obj* if present (works for lbftd ThermodynamicStates)."""
    if hasattr(obj, attr):
        object.__delattr__(obj, attr)
