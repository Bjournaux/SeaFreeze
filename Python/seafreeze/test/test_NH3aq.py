"""'NH3aq': the experimental aqueous-ammonia material.

The spline is built (SourceData/NH3aq) from a Redlich-Kister model fitted to ice
liquidus data.  These tests check, independently of the build code, that

  * the chemical potential of water carried by the spline is the fitted model
  * the pure-water limit is water1
  * melting curves and whichphase follow from it, with ice II on the liquidus
  * the calibration data are reproduced
  * uncalibrated properties are flagged, calibrated ones are not

The closed form used as reference is re-implemented here from the parameters
stored in the spline file, so a wrong G(P,T,m) construction cannot hide.
"""
import os
import pathlib
import warnings

import numpy as np
import pytest
from scipy.optimize import brentq

import seafreeze as sf
import seafreeze.seafreeze as sfm

R = 8.314462618
MW = sfm.mH2O_kgmol
NW = 1.0 / MW
M_NH3 = sfm.mNH3_kgmol
ICES = ('Ih', 'II', 'III', 'V', 'VI')
REPO = pathlib.Path(sfm.__file__).resolve().parents[2]
DATA = REPO / 'SourceData' / 'NH3aq' / 'liquidus_data.csv'


def wt_to_m(w):
    """NH3 mass percent -> molality [mol/kg]."""
    return (w / (100.0 - w)) / M_NH3


@pytest.fixture(scope='module')
def meta():
    sp = sfm._load_spline(sfm.defpath, 'NH3aq')
    rk = dict(zip(str(sp['rk_names']).split(), np.asarray(sp['rk_values'], float).ravel()))
    return sp, rk, float(sp['T_ref']), float(sp['P_ref']), float(sp['P_scale'])


def _pt(P, T):
    a = np.empty((1,), dtype=object); a[0] = (float(P), float(T))
    return a


def water_mu(P, T, T_ref=240.0):
    """Liquid-water reference [J/mol]: water1 above T_ref, constant-Cp continuation below."""
    if T >= T_ref:
        return sf.getProp(_pt(P, T), 'water1', sfm.defpath, 'G').G[0] * MW
    o = sf.getProp(_pt(P, T_ref), 'water1', sfm.defpath, 'G', 'S', 'Cp')
    g, s, cp = o.G[0] * MW, o.S[0] * MW, o.Cp[0] * MW
    return g - s * (T - T_ref) + cp * ((T - T_ref) - T * np.log(T / T_ref))


def muw_model(P, T, m, meta):
    """Chemical potential of water in the fitted Redlich-Kister model [J/mol]."""
    _, rk, T_ref, P_ref, P_scale = meta
    x2 = m / (NW + m); x1 = 1.0 - x2
    L = [rk[f'a{k}'] + rk[f'b{k}'] * (T_ref / T - 1.0) + rk[f'v{k}'] * (P - P_ref) / P_scale for k in range(3)]
    lng = x2**2 * (L[0] + L[1] * (4 * x1 - 1) + L[2] * (12 * x1**2 - 8 * x1 + 1))
    return water_mu(P, T, T_ref) + R * T * (np.log(x1) + lng)


def liquidus_model(P, m, meta, ices=ICES):
    """Stable liquidus of the model: (T, ice) with the highest melting temperature."""
    best = (np.nan, None)
    for ice in ices:
        def f(T):
            g = sf.getProp(_pt(P, T), ice, sfm.defpath, 'G').G[0] * MW
            return g - muw_model(P, T, m, meta)
        Tg = np.arange(201.0, 305.0, 2.0)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            r = np.array([f(t) for t in Tg])
        ok = np.isfinite(r[:-1]) & np.isfinite(r[1:]) & (np.signbit(r[:-1]) != np.signbit(r[1:]))
        if ok.any():
            i = np.where(ok)[0][-1]
            T = brentq(f, Tg[i], Tg[i + 1])
            if not np.isfinite(best[0]) or T > best[0]:
                best = (T, ice)
    return best


def _muw(P, T, m):
    a = np.empty((1,), dtype=object); a[0] = (float(P), float(T), float(m))
    return float(sf.getProp(a, 'NH3aq', sfm.defpath, 'muw').muw[0])


# ---------------------------------------------------------------------------
# Registration, domain, metadata
# ---------------------------------------------------------------------------
def test_registration_and_domain(meta):
    pd = sf.phases['NH3aq']
    assert (pd.solute, pd.nu, pd.experimental, pd.stitch) == ('NH3', 1, True, None)
    assert pd.MW == pytest.approx(17.03052e-3)
    rng = sf.phase_range('NH3aq')
    assert rng.P == (0.0, 1000.0) and rng.T == (200.0, 310.0) and rng.m == (0.0, 18.0)
    assert sf.solute_phase('NH3') == 'NH3aq'


def test_spline_file_documents_itself(meta):
    sp, rk, T_ref, P_ref, P_scale = meta
    assert (T_ref, P_ref, P_scale) == (240.0, 0.1, 1000.0)
    assert set(rk) == {f'{c}{k}' for c in 'abv' for k in range(3)}
    assert rk['a0'] < 0 and rk['b0'] < 0            # negative deviation, exothermic mixing
    for field in ('model', 'calibrated', 'not_calibrated', 'water_reference', 'nh3_reference', 'version'):
        assert len(str(sp[field])) > 5, field
    assert float(sp['MW'][1]) == pytest.approx(M_NH3) and int(sp['nu']) == 1
    assert 1.0 < float(sp['fit_rms_K']) < 1.5 and int(sp['fit_n_points']) == 104


# ---------------------------------------------------------------------------
# The calibrated quantity: chemical potential of water
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('P', [0.1, 55.0, 207.0, 433.0, 641.0, 977.0])
def test_muw_is_the_fitted_model(meta, P):
    """Scatter points between the building sites, on both sides of 240 K."""
    worst = 0.0
    for T in (203.3, 221.7, 239.0, 241.1, 263.9, 288.4, 307.7):
        for m in (1e-4, 0.013, 0.37, 2.9, 7.3, 12.6, 17.4):
            worst = max(worst, abs(_muw(P, T, m) - muw_model(P, T, m, meta)))
    assert worst < 0.1, f'P={P}: muw differs from the model by {worst:.3f} J/mol'    # ~0.005 K


def test_pure_water_limit_is_water1():
    P, T = np.array([0.1, 200.0, 600.0, 950.0]), np.array([240.0, 260.0, 300.0])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        s = sf.getProp(np.array([P, T, np.array([0.0])], dtype=object), 'NH3aq', sfm.defpath, 'G', 'rho', 'Cp', 'muw')
        w = sf.getProp(np.array([P, T], dtype=object), 'water1', sfm.defpath, 'G', 'rho', 'Cp')
    # 0.1 MPa is not a building site (0 and 10 MPa are): agreement is to interpolation accuracy
    np.testing.assert_allclose(s.G[:, :, 0], w.G, rtol=0, atol=1e-2)            # J/kg
    np.testing.assert_allclose(s.muw[:, :, 0], w.G * MW, rtol=0, atol=1e-3)     # J/mol
    np.testing.assert_allclose(s.rho[:, :, 0], w.rho, rtol=1e-5)
    # 240 K is where the water reference changes from water1 to its continuation: Cp is
    # continuous there but its slope is not, which the spline smooths slightly
    np.testing.assert_allclose(s.Cp[:, 1:, 0], w.Cp[:, 1:], rtol=1e-4)
    np.testing.assert_allclose(s.Cp[:, 0, 0], w.Cp[:, 0], rtol=2e-2)


def test_water_activity_matches_model(meta):
    P, T = 150.0, 255.0
    m = np.array([0.5, 3.0, 9.0, 16.0])
    out = sf.getProp(np.array([np.array([P]), np.array([T]), m], dtype=object), 'NH3aq', sfm.defpath, 'aw', 'phi')
    aw_model = np.array([np.exp((muw_model(P, T, mi, meta) - water_mu(P, T)) / (R * T)) for mi in m])
    np.testing.assert_allclose(np.ravel(out.aw), aw_model, rtol=2e-5)
    assert np.all(np.diff(np.ravel(out.aw)) < 0) and np.ravel(out.aw)[-1] < 0.75
    # strong negative deviation from ideality at high concentration: aw below the mole fraction
    assert np.ravel(out.aw)[-1] < NW / (NW + m[-1])


# ---------------------------------------------------------------------------
# Phase equilibria
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('ice, wt, P', [
    ('Ih', 5.0, 50.0), ('Ih', 11.0, 150.0), ('Ih', 20.0, 100.0),
    ('II', 11.0, 300.0), ('II', 20.0, 400.0), ('III', 5.0, 230.0),
    ('V', 11.0, 500.0), ('V', 20.0, 600.0), ('VI', 11.0, 800.0), ('VI', 20.0, 900.0),
])
def test_melting_curve_matches_model(meta, ice, wt, P):
    m = wt_to_m(wt)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        line = sf.phase_lines(ice, 'NH3aq', m=m)
    o = np.argsort(line.P)
    T_line = np.interp(P, line.P[o], line.T[o])
    T_model, _ = liquidus_model(P, m, meta, ices=(ice,))
    assert abs(T_line - T_model) < 0.05, f'{ice} {wt} wt% {P} MPa: {T_line:.3f} vs {T_model:.3f} K'


@pytest.mark.parametrize('wt, P, expected', [
    (8.0, 100.0, 'Ih'), (8.0, 230.0, 'III'), (8.0, 300.0, 'II'), (8.0, 500.0, 'V'), (8.0, 800.0, 'VI'),
    (11.0, 100.0, 'Ih'), (11.0, 230.0, 'II'), (11.0, 300.0, 'II'), (11.0, 500.0, 'V'), (11.0, 800.0, 'VI'),
    (20.0, 100.0, 'Ih'), (20.0, 300.0, 'II'), (20.0, 450.0, 'II'), (20.0, 600.0, 'V'), (20.0, 800.0, 'VI'),
])
def test_stable_ice_on_the_liquidus(meta, wt, P, expected):
    """Ice II replaces ice III on the liquidus as ammonia lowers the temperature."""
    m = wt_to_m(wt)
    T_liq, ice = liquidus_model(P, m, meta)
    assert ice == expected
    pts = np.empty((2,), dtype=object)
    pts[0], pts[1] = (P, T_liq - 0.5, m), (P, T_liq + 0.5, m)
    which = sf.whichphase(pts, 'NH3aq')
    assert sf.phasenum2phase(int(which[0])) == expected, f'{wt} wt% {P} MPa: below the liquidus'
    assert which[1] == 0, f'{wt} wt% {P} MPa: above the liquidus should be liquid'


def test_liquidus_depression_increases_with_ammonia():
    P = 100.0
    T = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for wt in (2.0, 8.0, 15.0, 22.0):
            line = sf.phase_lines('Ih', 'NH3aq', m=wt_to_m(wt))
            o = np.argsort(line.P)
            T.append(np.interp(P, line.P[o], line.T[o]))
        pure = sf.phase_lines('Ih', 'water1')
    op = np.argsort(pure.P)
    assert np.all(np.diff(T) < 0) and T[0] < np.interp(P, pure.P[op], pure.T[op])
    assert 40.0 < np.interp(P, pure.P[op], pure.T[op]) - T[-1] < 60.0      # ~49 K at 22 wt%


@pytest.mark.skipif(not DATA.exists(), reason='SourceData/NH3aq/liquidus_data.csv not available')
def test_reproduces_calibration_data():
    """The 104 liquidus points used for the fit, through getProp only."""
    import csv
    with open(DATA, newline='') as fh:
        d = list(csv.DictReader(line for line in fh if not line.startswith('#')))
    assert len(d) == 104
    dT = []
    for row in d:
        x = float(row['x_NH3']); m = x / (1 - x) * NW; P = float(row['P_MPa']); T = float(row['T_K'])
        ice = str(row['ice_fitted'])

        def f(t):
            return sf.getProp(_pt(P, t), ice, sfm.defpath, 'G').G[0] * MW - _muw(P, t, m)
        Tg = np.arange(202.0, 306.0, 4.0)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            r = np.array([f(t) for t in Tg])
        i = np.where(np.isfinite(r[:-1]) & np.isfinite(r[1:]) & (np.signbit(r[:-1]) != np.signbit(r[1:])))[0][-1]
        dT.append(brentq(f, Tg[i], Tg[i + 1]) - T)
    dT = np.array(dT)
    rms = float(np.sqrt(np.mean(dT**2)))
    assert 1.15 < rms < 1.30, f'RMS misfit {rms:.3f} K (fit: 1.20 K)'
    # the robust fit does not force a zero mean: the model runs about 0.4 K cold overall
    assert abs(dT.mean()) < 0.6 and np.abs(dT).max() < 5.0


# ---------------------------------------------------------------------------
# Experimental flag
# ---------------------------------------------------------------------------
PTM = np.array([np.array([100.0, 400.0]), np.array([250.0, 280.0]), np.array([1.0, 8.0])], dtype=object)


def _caught(fn):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        fn()
    return [x for x in w if issubclass(x.category, sf.ExperimentalMaterialWarning)]


def test_all_properties_warns_once_and_returns_full_set():
    w = _caught(lambda: sf.getProp(PTM, 'NH3aq'))
    assert len(w) == 1 and 'experimental' in str(w[0].message) and 'aw, muw, phi' in str(w[0].message)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        nh3, nacl = sf.getProp(PTM, 'NH3aq'), sf.getProp(PTM, 'NaClaq')
    assert set(vars(nh3)) == set(vars(nacl))            # same property set as NaCl
    for k, v in vars(nh3).items():
        assert np.all(np.isfinite(np.asarray(v, float))), k


def test_calibrated_properties_do_not_warn():
    assert _caught(lambda: sf.getProp(PTM, 'NH3aq', sfm.defpath, 'muw', 'aw', 'phi', 'm', 'xs', 'xw', 'f', 'P', 'T')) == []
    assert _caught(lambda: sf.phase_lines('Ih', 'NH3aq', m=5.0)) == []
    one_m = np.array([PTM[0], PTM[1], np.array([8.0])], dtype=object)     # whichphase grids take one molality
    assert _caught(lambda: sf.whichphase(one_m, 'NH3aq')) == []
    assert _caught(lambda: sf.getProp(PTM, 'NaClaq')) == []


def test_uncalibrated_properties_are_named():
    w = _caught(lambda: sf.getProp(PTM, 'NH3aq', sfm.defpath, 'rho', 'muw', 'Cp'))
    assert len(w) == 1 and 'Cp, rho' in str(w[0].message)
    assert sf.uncalibrated_properties('NH3aq', ('rho', 'muw', 'Cp', 'm')) == ['Cp', 'rho']
    assert sf.uncalibrated_properties('NH3aq', ('muw', 'aw')) == []
    assert sf.uncalibrated_properties('NaClaq', ('rho',)) == []
    assert len(sf.uncalibrated_properties('NH3aq')) == 1


def test_rho2P_works_but_warns():
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        rho = sf.getProp(np.array([np.array([300.0]), np.array([260.0]), np.array([5.0])], dtype=object),
                         'NH3aq', sfm.defpath, 'rho').rho.ravel()
    w = _caught(lambda: np.testing.assert_allclose(sf.rho2P(rho, 260.0, 'NH3aq', m=5.0), [300.0], atol=0.05))
    assert len(w) == 1 and 'rho' in str(w[0].message)


def test_wpd_overlay_runs():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        fig = sf.wpd(solute='NH3', m=[wt_to_m(10.0)])
    labels = [l.get_label() for l in fig.axes[0].get_lines()]
    assert any(lab.startswith('Ih-NH3aq') for lab in labels)
    plt.close('all')
