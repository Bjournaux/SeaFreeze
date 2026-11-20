"""Spline builder: function G(P, T[, m]) -> SeaFreeze Gibbs spline -> working material.

Four levels of evidence, from exact to end-to-end:
  1. polynomials inside the spline space are reproduced to round-off
  2. an existing SeaFreeze solution spline (NaClaq_LP) is rebuilt coefficient for
     coefficient, saved, reloaded, registered and gives the same getProp output
  3. a solution-like function with the singular m ln m mixing term reaches a
     stated accuracy on the water chemical potential, and improves on refinement
  4. a synthetic ideal binary built from a function is registered and yields the
     thermodynamically expected water activity and ice Ih melting-point depression
"""
import os
import warnings

import numpy as np
import pytest
from scipy.optimize import brentq
from scipy.special import xlogy

import seafreeze as sf
import seafreeze.seafreeze as sfm
from seafreeze import splinebuilder as sb
from mlbspline import load
from mlbspline.eval import evalMultivarSpline

R = 8.314462618
MWV = sfm.mH2O_kgmol
NW = 1.0 / MWV


def _eval(sp, axes, der):
    x = np.empty(len(axes), dtype=object)
    for i, a in enumerate(axes):
        x[i] = np.asarray(a, float)
    return evalMultivarSpline(sp, x, list(der), allowExtrapolations=True)


def _water_G(P, T):
    """water1 specific Gibbs energy on a (nP, nT) grid."""
    return sf.getProp(np.array([np.asarray(P, float), np.asarray(T, float)], dtype=object),
                      'water1', sfm.defpath, 'G').G


def _ideal_binary(P_sites, T_sites, MWs):
    """G(P, T, m) of an ideal solution of a non-electrolyte in water1 [J/kg solution]."""
    Gw = _water_G(P_sites, T_sites)

    def G(P, T, m):
        mus0 = 1.0e3 + 20.0 * P - 30.0 * (T - 273.0)              # arbitrary smooth solute term
        mix = R * T * (xlogy(m, m / (m + NW)) + NW * np.log(NW / (m + NW)))
        return (Gw[:, :, None] + m * mus0 + mix) / (1.0 + MWs * m)
    return G


# ---------------------------------------------------------------------------
# 1. Exactness on polynomials
# ---------------------------------------------------------------------------
def test_polynomial_reproduced_exactly_3d():
    def G(P, T, m):
        p, t = P / 500.0, (T - 300.0) / 50.0
        return 1e4 * (1 + p - 0.3 * p**2 + 0.05 * p**5) * (2 - t + 0.4 * t**3 - 0.1 * t**5) * (1 + m - 0.2 * m**3)
    sp = sb.build_gibbs_spline(G, np.linspace(0, 1000, 9), np.linspace(250, 350, 8), np.linspace(0, 6, 7))
    assert sp['form'] == 'B-' and tuple(sp['order']) == (6, 6, 4)
    assert sp['coefs'].shape == tuple(sp['number']) == (9, 8, 7)
    Pt, Tt, mt = np.array([13.0, 487.0, 951.0]), np.array([251.5, 301.0, 349.0]), np.array([0.01, 2.2, 5.9])
    mesh = (Pt[:, None, None], Tt[None, :, None], mt[None, None, :])
    np.testing.assert_allclose(_eval(sp, (Pt, Tt, mt), (0, 0, 0)), G(*mesh), rtol=1e-11)
    h = 1e-3
    dP = (G(mesh[0] + h, mesh[1], mesh[2]) - G(mesh[0] - h, mesh[1], mesh[2])) / (2 * h)
    dm = (G(mesh[0], mesh[1], mesh[2] + h) - G(mesh[0], mesh[1], mesh[2] - h)) / (2 * h)
    np.testing.assert_allclose(_eval(sp, (Pt, Tt, mt), (1, 0, 0)), dP, rtol=1e-6)
    np.testing.assert_allclose(_eval(sp, (Pt, Tt, mt), (0, 0, 1)), dm, rtol=1e-6)


def test_polynomial_reproduced_exactly_2d_pure_substance():
    G = lambda P, T: 5e3 + 900.0 * P - 0.2 * P**2 - 3.0 * T * np.log(300.0) + 1e-3 * P * T - 0.02 * T**2
    sp = sb.build_gibbs_spline(G, np.linspace(0, 500, 8), np.linspace(200, 300, 9))
    assert len(sp['knots']) == 2 and tuple(sp['order']) == (6, 6)
    Pt, Tt = np.array([3.0, 222.0, 499.0]), np.array([201.0, 250.5, 299.0])
    np.testing.assert_allclose(_eval(sp, (Pt, Tt), (0, 0)), G(Pt[:, None], Tt[None, :]), rtol=1e-12)
    err = sb.check_gibbs_spline(sp, G, Pt, Tt)
    assert err['G'] < 1e-6 and err['Cp'] < 1e-3 and err['rho'] < 1e-3


# ---------------------------------------------------------------------------
# 2. Round trip of a real SeaFreeze solution spline
# ---------------------------------------------------------------------------
@pytest.fixture(scope='module')
def nacl_roundtrip(tmp_path_factory):
    sp0 = sfm._load_spline(sfm.defpath, 'NaClaq_LP')
    sites = [sb.greville_sites(k, o) for k, o in zip(sp0['knots'], sp0['order'])]
    values = _eval(sp0, sites, (0, 0, 0))
    sp1 = sb.build_gibbs_spline(values, *sites, order=sp0['order'], knots=list(sp0['knots']))
    d = tmp_path_factory.mktemp('splines')
    fn = sb.save_gibbs_spline(sp1, os.path.join(d, 'rt', 'rt.mat'), MW=float(sp0['MW'][1]), nu=2,
                              model='round trip of NaClaq_LP', calibrated='all')
    return sp0, sp1, str(d), fn


def test_roundtrip_coefficients(nacl_roundtrip):
    sp0, sp1, _, _ = nacl_roundtrip
    for a, b in zip(sp0['knots'], sp1['knots']):
        np.testing.assert_array_equal(a, b)
    scale = np.abs(sp0['coefs']).max()
    assert np.abs(sp1['coefs'] - sp0['coefs']).max() < 1e-11 * scale


def test_roundtrip_file_is_loadable_and_keeps_metadata(nacl_roundtrip):
    sp0, _, _, fn = nacl_roundtrip
    sp2 = load.loadSpline(fn, 'sp')
    assert sp2['form'] == 'B-' and int(sp2['dim']) == 1
    np.testing.assert_array_equal(sp2['number'], sp0['number'])
    np.testing.assert_array_equal(sp2['order'], sp0['order'])
    for a, b in zip(sp0['knots'], sp2['knots']):
        np.testing.assert_array_equal(a, b)
    assert float(sp2['MW']) == float(sp0['MW'][1]) and float(sp2['nu']) == 2.0
    assert str(sp2['model']) == 'round trip of NaClaq_LP' and str(sp2['calibrated']) == 'all'


def test_roundtrip_material_matches_getProp(nacl_roundtrip):
    sp0, _, d, _ = nacl_roundtrip
    # registry molar mass as for the built-in NaCl materials (58.44e-3; the spline file carries 0.058443)
    sf.register_solution('TESTrt', 'TESTRT', sfm.mNaCl_kgmol, 2, 0.0002, spline=('rt', 'rt.mat'))
    try:
        grid = np.array([np.array([0.1, 200.0, 600.0, 900.0]), np.array([240.0, 280.0, 350.0, 450.0]),
                         np.array([0.0, 0.5, 3.0, 6.0])], dtype=object)
        scat = np.empty((3,), dtype=object)
        scat[0], scat[1], scat[2] = (50.0, 300.0, 1.0), (450.0, 260.0, 4.0), (950.0, 400.0, 0.1)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            for PTm in (grid, scat):
                a, b = sf.getProp(PTm, 'TESTrt', d), sf.getProp(PTm, 'NaClaq_LP')
                assert set(vars(a)) == set(vars(b))
                for k, v in vars(a).items():
                    A, B = np.asarray(v, float), np.asarray(getattr(b, k), float)
                    np.testing.assert_array_equal(np.isfinite(A), np.isfinite(B), err_msg=k)
                    np.testing.assert_allclose(A, B, rtol=1e-6, atol=1e-9, equal_nan=True, err_msg=k)
            assert sf.phase_range('TESTrt', d) == sf.phase_range('NaClaq_LP')
    finally:
        sf.unregister_solution('TESTrt')


# ---------------------------------------------------------------------------
# 3. Dilute limit: the m ln m term
# ---------------------------------------------------------------------------
def _muw_error(n_m, power):
    Pg, Tg = np.linspace(0.1, 300.0, 12), np.linspace(245.0, 300.0, 12)
    sp = sb.build_gibbs_spline(_ideal_binary(Pg, Tg, 0.017), Pg, Tg, sb.dilute_sites(6.0, n_m, power))
    mt = np.array([1e-4, 1e-3, 0.01, 0.05, 0.2, 0.7, 1.3, 2.9, 4.4, 5.9])       # off-site molalities
    G0, G1 = _eval(sp, (Pg, Tg, mt), (0, 0, 0)), _eval(sp, (Pg, Tg, mt), (0, 0, 1))
    muw = MWV * (G0 - (1 + 0.017 * mt) * mt * G1)
    exact = _water_G(Pg, Tg)[:, :, None] * MWV + R * Tg[None, :, None] * np.log(NW / (mt + NW))
    return float(np.abs(muw - exact).max())


def test_dilute_sites_shape():
    s = sb.dilute_sites(6.0, 11, 3.0)
    assert s[0] == 0.0 and s[-1] == 6.0 and np.all(np.diff(s) > 0)
    assert s[1] - s[0] < 0.01 * (s[-1] - s[-2])            # much denser at the dilute end


def test_water_potential_accuracy_and_convergence():
    """muw of an ideal solution, exact value G_w + RT ln(x_w): J/mol of H2O.
    0.05 J/mol is about 0.002 K on a melting curve."""
    coarse, fine = _muw_error(20, 3.0), _muw_error(40, 3.0)
    assert fine < 0.05, f'muw error {fine:.4f} J/mol'
    assert fine < 0.5 * coarse, f'no convergence: {coarse:.4f} -> {fine:.4f} J/mol'
    assert _muw_error(40, 3.0) < _muw_error(40, 1.0), 'clustering near m = 0 should help'


# ---------------------------------------------------------------------------
# 4. End to end: function -> spline -> registered binary -> melting curve
# ---------------------------------------------------------------------------
@pytest.fixture(scope='module')
def ideal_material(tmp_path_factory):
    MWs = 0.017
    Pg, Tg, mg = np.linspace(0.1, 300.0, 25), np.linspace(240.0, 300.0, 25), sb.dilute_sites(6.0, 40, 3.0)
    G = _ideal_binary(Pg, Tg, MWs)
    sp = sb.build_gibbs_spline(G, Pg, Tg, mg)
    d = tmp_path_factory.mktemp('splines_ideal')
    sb.save_gibbs_spline(sp, os.path.join(d, 'ideal', 'ideal.mat'), MW=MWs, nu=1, model='ideal binary')
    # the ices and water1 are read from the same splines directory
    for sub in ('ice_Ih', 'ice_II', 'ice_III', 'ice_V', 'ice_VI', 'ice_VII_X_French', 'water_Bollengier'):
        os.symlink(os.path.join(sfm.defpath, sub), os.path.join(d, sub))
    sf.register_solution('TESTideal', 'TESTID', MWs, 1, 0.0002, spline=('ideal', 'ideal.mat'),
                         experimental=True)
    yield str(d), (Pg, Tg, mg), G, sp, MWs
    sf.unregister_solution('TESTideal')


def test_check_gibbs_spline_reports_small_errors(ideal_material):
    _, (Pg, Tg, _), G, sp, MWs = ideal_material
    err = sb.check_gibbs_spline(sp, G, Pg, Tg, m=np.array([0.003, 0.4, 2.7, 5.5]), MW=MWs)
    assert set(err) == {'G', 'rho', 'Cp', 'muw'}
    assert err['muw'] < 0.05      # J/mol of H2O
    assert err['G'] < 0.5         # J/kg
    # 'rho' and 'Cp' are not asserted here: this G takes water from a table on the
    # building sites, so it cannot be differentiated in P or T by finite differences.
    # They are checked on analytic functions in the polynomial tests.


def test_registered_binary_water_activity_is_ideal(ideal_material):
    """For a non-electrolyte (nu = 1) aw must equal x_w.  This caught lbftd
    hard-coding nu = 2 in the water activity (it returned x_w squared)."""
    d, _, _, _, _ = ideal_material
    m = np.array([0.0, 0.5, 2.0, 5.0])
    out = sf.getProp(np.array([np.array([0.1, 150.0]), np.array([260.0, 290.0]), m], dtype=object),
                     'TESTideal', d, 'aw', 'muw', 'xw', 'm', 'f')
    # aw is referenced to the solution at the dilute cutoff (0.0002 mol/kg), not to m = 0,
    # which biases ln(aw) by about 5e-5 here; muw has no such bias and is what phase
    # equilibria use.
    np.testing.assert_allclose(out.aw, np.broadcast_to(NW / (m + NW), out.aw.shape), rtol=1e-4)
    Tg = np.array([260.0, 290.0])
    muw_exact = (_water_G([0.1, 150.0], Tg)[:, :, None] * MWV
                 + R * Tg[None, :, None] * np.log(NW / (m + NW)))
    assert np.abs(out.muw - muw_exact).max() < 0.05          # J/mol of H2O
    np.testing.assert_allclose(out.xw, np.broadcast_to(NW / (m + NW), out.xw.shape), rtol=1e-12)   # nu = 1
    np.testing.assert_allclose(out.f, np.broadcast_to(1 + 0.017 * m, out.f.shape), rtol=1e-12)
    assert sf.phases['TESTideal'].experimental and sf.is_solution('TESTideal')


def test_registered_binary_melting_depression(ideal_material):
    """Ice Ih liquidus of the ideal binary = solution of G_Ih = G_w + RT ln(x_w)."""
    d, _, _, _, _ = ideal_material
    m_val = 2.0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        line = sf.phase_lines('Ih', 'TESTideal', m=m_val, path=d)
        pure = sf.phase_lines('Ih', 'water1')
    order = np.argsort(line.P)

    def T_expected(P):
        def f(T):
            pt = np.empty((1,), dtype=object); pt[0] = (P, T)
            gi = sf.getProp(pt, 'Ih', sfm.defpath, 'G').G[0] * MWV
            gw = sf.getProp(pt, 'water1', sfm.defpath, 'G').G[0] * MWV
            return gi - gw - R * T * np.log(NW / (m_val + NW))
        return brentq(f, 245.0, 275.0)

    for P in (20.0, 100.0, 180.0):
        T_line = np.interp(P, line.P[order], line.T[order])
        assert abs(T_line - T_expected(P)) < 0.05, f'P={P}: {T_line:.3f} vs {T_expected(P):.3f}'
    # and it is a depression of the right size (ideal cryoscopic: ~1.86 K per mol/kg at 1 bar)
    op = np.argsort(pure.P)
    dT = np.interp(20.0, pure.P[op], pure.T[op]) - np.interp(20.0, line.P[order], line.T[order])
    assert 3.0 < dT < 4.5
    # whichphase agrees with the melting line on both sides of it
    T_m = np.interp(100.0, line.P[order], line.T[order])
    pts = np.empty((2,), dtype=object)
    pts[0], pts[1] = (100.0, T_m - 2.0, m_val), (100.0, T_m + 2.0, m_val)
    np.testing.assert_array_equal(sf.whichphase(pts, 'TESTideal', d), [1, 0])


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------
def test_build_validation():
    G = lambda P, T, m: P + T + m
    P, T, m = np.linspace(0, 1, 8), np.linspace(1, 2, 8), np.linspace(0, 1, 6)
    with pytest.raises(ValueError, match='at least 6 sites'):
        sb.build_gibbs_spline(G, P[:4], T, m)
    with pytest.raises(ValueError, match='strictly increasing'):
        sb.build_gibbs_spline(G, P[::-1], T, m)
    with pytest.raises(ValueError, match='order must have 3'):
        sb.build_gibbs_spline(G, P, T, m, order=(6, 6))
    with pytest.raises(ValueError, match='not finite'):
        sb.build_gibbs_spline(lambda P, T, m: np.log(m) + P + T, P, T, m)
    with pytest.raises(ValueError, match='expected'):
        sb.build_gibbs_spline(np.zeros((3, 3, 3)), P, T, m)


def test_save_rejects_clashing_metadata(tmp_path):
    sp = sb.build_gibbs_spline(lambda P, T: P + T, np.linspace(0, 1, 7), np.linspace(1, 2, 7))
    with pytest.raises(ValueError, match='clashes'):
        sb.save_gibbs_spline(sp, str(tmp_path / 'x.mat'), coefs='oops')
