"""Tests for the Helmholtz liquid 'water3' (psi surface) and lbftd.evalHelmholtz.

References (Matlab/test/fixtures/, written by Matlab/test/gen_water3_reference.m):
  psi_reference.mat             psiH2O_val.m (lbf-thermo) at (rho,T) states
  water3_getprop_reference.mat  Matlab SF_getprop(...,'water3') on a (P,T) grid
  water_F_test.mat              plain F(rho,T) spline (IAPWS-95 fit), generic path
"""
import os
import warnings

import numpy as np
import pytest
from scipy.io import loadmat

import seafreeze as sf
from seafreeze.seafreeze import defpath
from mlbspline import load
from lbftd import evalHelmholtz as eh

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.normpath(os.path.join(HERE, '..', '..', '..', 'Matlab', 'test', 'fixtures'))


def _fix(name):
    p = os.path.join(FIX, name)
    if not os.path.exists(p):
        pytest.skip(f'{name} missing; run Matlab/test/gen_water3_reference.m')
    return loadmat(p, squeeze_me=True, struct_as_record=False)


def _scatter(a, b):
    out = np.empty(len(a), dtype=object)
    for i, (x, y) in enumerate(zip(a, b)):
        out[i] = (float(x), float(y))
    return out


def _grid(a, b):
    return np.array([np.asarray(a, float), np.asarray(b, float)], dtype=object)


def _relerr(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    return np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-300))


@pytest.fixture(scope='module')
def sp3():
    return load.loadSpline(os.path.join(defpath, 'water_psi2026', 'water_psi2026.mat'), 'sp')


# --------------------------------------------------------------------------
# 1. evaluator vs psiH2O_val
# --------------------------------------------------------------------------
PAIRS = [('P', 'P'), ('G', 'G'), ('S', 'S'), ('U', 'U'), ('H', 'H'), ('A', 'F'), ('Cp', 'cp'),
         ('Cv', 'cv'), ('Kt', 'Kt'), ('Ks', 'Ks'), ('alpha', 'alpha')]


def test_psi_scatter_vs_psiH2O_val(sp3):
    ref = _fix('psi_reference.mat')['ref']
    o = eh.evalHelmholtzScatter(sp3, _scatter(ref.rho, ref.T), rhoT=True)
    for a, b in PAIRS:
        x = getattr(o, a)
        y = np.real(getattr(ref.out, b))
        ok = np.isfinite(y)
        assert np.array_equal(np.isnan(x), np.isnan(y)), a
        assert _relerr(x[ok], y[ok]) < 1e-9, a
    w = ref.out.w
    okw = np.isfinite(w) & (np.imag(w) == 0)
    assert _relerr(o.vel[okw], np.real(w[okw])) < 1e-9
    assert np.all(np.isnan(o.vel[np.imag(w) != 0]))


def test_psi_grid_and_edges(sp3):
    ref = _fix('psi_reference.mat')['ref']
    g = eh.evalHelmholtzGrid(sp3, _grid(ref.grid.rho, ref.grid.T), 'P', 'G', rhoT=True)
    assert g.P.shape == (ref.grid.rho.size, ref.grid.T.size)
    assert _relerr(g.P, ref.grid.out.P) < 1e-9
    assert _relerr(g.G, ref.grid.out.G) < 1e-9
    e = eh.evalHelmholtzScatter(sp3, _scatter(ref.edge.rho, ref.edge.T), 'P', 'G', rhoT=True)
    assert _relerr(e.P[0], ref.edge.out.P[0]) < 1e-9          # below-floor virial continuation
    assert np.all(np.isnan(e.P[1:])) and np.all(np.isnan(e.G[1:]))


# --------------------------------------------------------------------------
# 2. getProp (P,T): Matlab parity, identities, round trip
# --------------------------------------------------------------------------
def test_getProp_water3_vs_matlab():
    w3 = _fix('water3_getprop_reference.mat')['w3']
    o = sf.getProp(_grid(w3.P, w3.T), 'water3')
    for k in ['rho', 'G', 'S', 'U', 'H', 'A', 'Cp', 'Cv', 'Kt', 'Kp', 'Ks', 'alpha', 'vel',
              'Js', 'gamma_Gruneisen']:
        x = np.asarray(getattr(o, k), float)
        y = np.asarray(getattr(w3.grid, k), float)
        assert np.array_equal(np.isnan(x), np.isnan(y)), k
        ok = np.isfinite(y)
        assert _relerr(x[ok], y[ok]) < (1e-4 if k == 'Kp' else 1e-8), k


def test_getProp_water3_grid_scatter_and_identities():
    P = np.array([0.1, 50, 500, 1500.]); T = np.array([260., 300, 350, 400])
    g = sf.getProp(_grid(P, T), 'water3')
    assert g.rho.shape == (4, 4)
    Pm, Tm = np.meshgrid(P, T, indexing='ij')
    s = sf.getProp(_scatter(Pm.ravel(), Tm.ravel()), 'water3')
    assert _relerr(s.rho, g.rho.ravel()) < 1e-12
    assert _relerr(s.Cp - s.Cv, s.T * s.alpha ** 2 * s.Kt * 1e6 / s.rho) < 1e-8
    assert _relerr(s.Ks / s.Kt, s.Cp / s.Cv) < 1e-10
    assert _relerr(s.vel ** 2, s.Ks * 1e6 / s.rho) < 1e-10
    assert _relerr(s.G, s.U - s.T * s.S + s.P * 1e6 / s.rho) < 1e-8
    # round trip through rhoT
    b = sf.getProp(_scatter(s.rho, s.T), 'water3', defpath, 'P', rhoT=True)
    assert np.max(np.abs(b.P - Pm.ravel())) < 1e-7


def test_getProp_water3_ambient_and_domain():
    a = sf.getProp(_scatter([0.101325], [298.15]), 'water3')
    assert abs(a.rho[0] - 997.05) < 0.05
    assert abs(a.Cp[0] - 4181.5) < 5
    assert abs(a.vel[0] - 1496.7) < 1
    z = sf.getProp(_scatter([100, 1e7], [200, 300]), 'water3', defpath, 'rho', 'G')
    assert np.all(np.isnan(z.rho)) and np.all(np.isnan(z.G))


def test_branch_selection_water3():
    pts = _scatter([1e-3, 0.1, 0.1, 10], [300, 400, 300, 400])
    v = sf.getProp(pts, 'water3', defpath, 'rho', 'G')
    assert v.rho[0] < 0.01 and v.rho[1] < 1 and v.rho[2] > 990 and v.rho[3] > 930
    lq = sf.getProp(_scatter([1e-3, 0.1], [300, 400]), 'water3', defpath, 'rho', 'G', branch='liquid')
    assert np.all(lq.rho > 930) and np.all(lq.G > v.G[:2])
    vp = sf.getProp(_scatter([0.1], [300]), 'water3', defpath, 'rho', 'G', branch='vapor')
    assert vp.rho[0] < 1 and vp.G[0] > v.G[2]


def test_rhoT_refused_for_gibbs():
    with pytest.raises(ValueError):
        sf.getProp(_scatter([1000], [300]), 'water1', defpath, 'rho', rhoT=True)


def test_rho2P_and_range_water3():
    P = sf.rho2P([997.047, 1100.0, 1200.0], [298.15, 300.0, 350.0], 'water3')
    back = sf.getProp(_scatter(P, [298.15, 300.0, 350.0]), 'water3', defpath, 'rho')
    assert np.max(np.abs(back.rho - [997.047, 1100.0, 1200.0])) < 1e-6
    r = sf.phase_range('water3')
    assert r.rho is not None and r.T[0] < 240 and r.P[1] > 2300


# --------------------------------------------------------------------------
# 3. phase equilibria: reference state consistent with the ice splines
# --------------------------------------------------------------------------
def test_whichphase_water3_matches_water1():
    PT = _grid([0.1, 100, 300, 800, 1500], [250, 260, 270, 276, 300, 330])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        assert np.array_equal(sf.whichphase(PT, 'water3'), sf.whichphase(PT))


def test_melting_Ih_water3():
    r = sf.phase_lines('Ih', 'water3', P=np.array([0.1, 0.101325, 1.0]),
                       T=np.arange(270, 276, 0.005))
    T01 = np.interp(0.101325, np.sort(r.P), r.T[np.argsort(r.P)])
    assert abs(T01 - 273.152) < 0.01


def test_melting_curves_water3_vs_water1():
    for ice, tol in [('Ih', 0.03), ('III', 0.03), ('V', 0.08), ('VI', 0.5)]:
        a = sf.phase_lines(ice, 'water3', segment='stable')
        b = sf.phase_lines(ice, 'water1', segment='stable')
        ia, ib = np.argsort(a.P), np.argsort(b.P)
        Pc = np.linspace(max(a.P.min(), b.P.min()), min(a.P.max(), b.P.max()), 30)
        dT = np.interp(Pc, a.P[ia], a.T[ia]) - np.interp(Pc, b.P[ib], b.T[ib])
        assert np.max(np.abs(dT)) < tol, ice


# --------------------------------------------------------------------------
# 4. generic F(rho,T) spline path
# --------------------------------------------------------------------------
def test_plain_F_rhoT_spline():
    p = os.path.join(FIX, 'water_F_test.mat')
    if not os.path.exists(p):
        pytest.skip('water_F_test.mat missing; run Matlab/test/gen_helmholtz_fixture.m')
    spF = load.loadSpline(p, 'sp')
    q = eh.evalHelmholtzGrid(spF, _grid([1, 100, 1000], [280, 300, 400]), 'rho', 'Cp')
    assert np.all(np.isfinite(q.rho))
    Pm, Tm = np.meshgrid([1, 100, 1000], [280, 300, 400], indexing='ij')
    b = eh.evalHelmholtzScatter(spF, _scatter(q.rho.ravel(), Tm.ravel()), 'P', rhoT=True)
    assert np.max(np.abs(b.P - Pm.ravel())) < 1e-7


# --------------------------------------------------------------------------
# 5. saturation / sublimation
# --------------------------------------------------------------------------
def test_saturation_vs_iapws_aux():
    from seafreeze import coexistence as cx
    T = np.array([273.16, 300, 373.124, 450, 550, 620, 640])
    s = sf.saturation(T)
    assert np.max(np.abs(s.P / cx.psat_iapws_aux(T) - 1)) < 3e-4
    assert abs(s.P[2] - 0.101325) < 2e-5
    assert np.all(s.rho_A > s.rho_B)


def test_saturation_smooth_near_Tc():
    T = np.arange(620, 646.6, 0.25)
    s = sf.saturation(T)
    assert np.all(np.isfinite(s.P))
    assert np.all(np.diff(s.P) > 0)                      # monotonic vapour pressure
    assert np.all(np.diff(s.rho_B) > 0)
    # stage5_18f: a small (dP/drho)_T loop straddles the liquid binodal at
    # ~636.5-637.5 K (psiEOS LIMITS: sign changes of (dP/drho)_T at 600-642 K),
    # so rho' steps by ~10 kg/m3 there.  Monotonic everywhere else.
    ok = (T[:-1] < 636.5) | (T[:-1] > 637.75)
    assert np.all(np.diff(s.rho_A)[ok] < 0)
    assert np.max(np.abs(np.diff(s.rho_A)[~ok])) < 15


def test_sublimation_vs_R1408_and_NIST():
    from seafreeze import coexistence as cx
    T = np.array([230, 240, 250, 260, 270, 273.16])
    s = sf.sublimation(T)
    assert np.max(np.abs(s.P / cx.psub_iapws(T) - 1)) < 2e-4
    # below 230 K: NaN by default, dilute-vapour extension on request
    assert np.all(np.isnan(sf.sublimation([175.0, 200.0]).P))
    from seafreeze.test.water3_vapor_figures import BIELSKA2013
    Tb, pb, ub = BIELSKA2013.T
    e = sf.sublimation(Tb, dilute_extension=True)
    # every NIST Bielska et al. (2013) point within 3 sigma
    assert np.all(np.abs(e.P * 1e6 - pb) < 3 * ub)


def test_coexistence_vs_matlab():
    w3 = _fix('water3_getprop_reference.mat')['w3']
    if not hasattr(w3, 'sat'):
        pytest.skip('reference predates SF_coexistence; rerun gen_water3_reference.m')
    s = sf.saturation(w3.sat.T)
    assert _relerr(s.P, w3.sat.P) < 1e-9
    b = sf.sublimation(w3.sub.T, dilute_extension=True)
    assert _relerr(b.P, w3.sub.P) < 1e-9
