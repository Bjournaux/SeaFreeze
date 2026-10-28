"""Domain-boundary and edge-case tests for getProp.

* Points outside the spline domain return NaN, never crash, and keep shape —
  including the all-points-out case that used to return an empty object.
* Mixed in/out scatter sets return NaN only where expected.
* The NaClaq LP+HP stitch is an exact cosine blend of the two sub-splines,
  and is continuous at both taper boundaries.
"""

import numpy as np
import pytest

import seafreeze as sf
from seafreeze.seafreeze import defpath as DEFPATH
from seafreeze.seafreeze import _STITCH_P_LO, _STITCH_P_HI

OUT_OF_DOMAIN = {
    # material: a (P, T[, m]) point that is outside its spline domain
    'Ih':             (3000.0, 250.0),
    'VI':             (5000.0, 250.0),
    'VII_X_French':   (300.0, 260.0),
    'water1':         (5000.0, 300.0),
    'NaClaq_HP':      (300.0, 260.0, 1.0),
    'NaClaq_LP':      (3000.0, 300.0, 1.0),
    'NaClaq':         (50000.0, 300.0, 1.0),
}
IN_DOMAIN = {
    'Ih':             (100.0, 250.0),
    'VI':             (1500.0, 280.0),
    'VII_X_French':   (5000.0, 300.0),
    'water1':         (500.0, 300.0),
    'NaClaq_HP':      (2000.0, 600.0, 1.0),
    'NaClaq_LP':      (300.0, 300.0, 1.0),
    'NaClaq':         (750.0, 300.0, 1.0),
}


def _scatter(*pts):
    a = np.empty(len(pts), object)
    for i, p in enumerate(pts):
        a[i] = tuple(float(x) for x in p)
    return a


@pytest.mark.parametrize('material', sorted(OUT_OF_DOMAIN))
def test_all_points_out_of_domain_all_mode(material):
    out = vars(sf.getProp(_scatter(OUT_OF_DOMAIN[material]), material))
    assert out, f'{material}: empty output for out-of-domain point'
    for k, v in out.items():
        v = np.asarray(v, float)
        assert v.shape == (1,), f'{material}/{k}: shape {v.shape}'
        if k in ('P', 'T', 'm', 'xs', 'xw', 'f', 'Vw'):
            continue  # echoes / m-only quantities / extrapolated Vw are defined
        assert np.isnan(v).all(), f'{material}/{k} = {v} should be NaN out of domain'


@pytest.mark.parametrize('material', sorted(OUT_OF_DOMAIN))
def test_all_points_out_of_domain_selective(material):
    out = vars(sf.getProp(_scatter(OUT_OF_DOMAIN[material]), material, DEFPATH, 'rho', 'Js'))
    assert set(out) == {'rho', 'Js'}
    assert np.isnan(out['rho']).all() and np.isnan(out['Js']).all()


@pytest.mark.parametrize('material', sorted(OUT_OF_DOMAIN))
def test_mixed_in_and_out_of_domain(material):
    out = sf.getProp(_scatter(OUT_OF_DOMAIN[material], IN_DOMAIN[material]), material,
                     DEFPATH, 'rho', 'G')
    for k in ('rho', 'G'):
        v = np.asarray(getattr(out, k), float)
        assert v.shape == (2,)
        assert np.isnan(v[0]), f'{material}/{k}[0] should be NaN'
        assert np.isfinite(v[1]), f'{material}/{k}[1] should be finite'


def test_grid_partially_out_of_domain_is_nan_only_there():
    PT = np.empty(2, object)
    PT[0] = np.array([100.0, 2000.0])     # 2000 MPa is outside ice Ih
    PT[1] = np.array([240.0, 260.0])
    rho = np.asarray(sf.getProp(PT, 'Ih', DEFPATH, 'rho').rho, float)
    assert rho.shape == (2, 2)
    assert np.isfinite(rho[0]).all()
    assert np.isnan(rho[1]).all()


# ---------------------------------------------------------------------------
# NaClaq stitching
# ---------------------------------------------------------------------------
_STITCH_PROPS = ('G', 'rho', 'Cp', 'alpha', 'Kt', 'muw', 'mus', 'Va', 'phi', 'aw')


def test_stitch_is_lp_only_below_taper():
    pt = _scatter((200.0, 300.0, 2.0))
    st = vars(sf.getProp(pt, 'NaClaq', DEFPATH, *_STITCH_PROPS))
    lp = vars(sf.getProp(pt, 'NaClaq_LP', DEFPATH, *_STITCH_PROPS))
    for k in _STITCH_PROPS:
        np.testing.assert_allclose(st[k], lp[k], rtol=1e-12, err_msg=k)


def test_stitch_is_hp_only_above_taper():
    pt = _scatter((3000.0, 300.0, 2.0))
    st = vars(sf.getProp(pt, 'NaClaq', DEFPATH, *_STITCH_PROPS))
    hp = vars(sf.getProp(pt, 'NaClaq_HP', DEFPATH, *_STITCH_PROPS))
    for k in _STITCH_PROPS:
        np.testing.assert_allclose(st[k], hp[k], rtol=1e-12, err_msg=k)


@pytest.mark.parametrize('P', [550.0, 750.0, 950.0])
def test_stitch_is_exact_cosine_blend_in_taper(P):
    """Inside [499.9, 1001] MPa every property is w*LP + (1-w)*HP with the
    documented cosine weight.  Guards both the weight formula and the
    index bookkeeping in _nacl_stitch."""
    w = 0.5 * (1.0 + np.cos(np.pi * (P - _STITCH_P_LO) / (_STITCH_P_HI - _STITCH_P_LO)))
    pt = _scatter((P, 300.0, 2.0))
    st = vars(sf.getProp(pt, 'NaClaq', DEFPATH, *_STITCH_PROPS))
    lp = vars(sf.getProp(pt, 'NaClaq_LP', DEFPATH, *_STITCH_PROPS))
    hp = vars(sf.getProp(pt, 'NaClaq_HP', DEFPATH, *_STITCH_PROPS))
    for k in _STITCH_PROPS:
        expected = w * np.asarray(lp[k]) + (1.0 - w) * np.asarray(hp[k])
        np.testing.assert_allclose(st[k], expected, rtol=1e-10, err_msg=f'{k} at P={P}')


@pytest.mark.parametrize('P_edge', [_STITCH_P_LO, _STITCH_P_HI])
def test_stitch_continuous_at_taper_boundaries(P_edge):
    eps = 1e-4
    pts = _scatter((P_edge - eps, 300.0, 2.0), (P_edge + eps, 300.0, 2.0))
    out = vars(sf.getProp(pts, 'NaClaq', DEFPATH, 'G', 'rho', 'muw'))
    for k, v in out.items():
        v = np.asarray(v, float)
        assert np.isfinite(v).all(), f'{k} not finite at taper edge {P_edge}'
        rel = abs(v[1] - v[0]) / max(abs(v[0]), 1e-300)
        assert rel < 1e-6, f'{k} jumps across P={P_edge}: {v}'


def test_stitch_grid_mixes_zones_correctly():
    """Grid with P axis spanning LP-only, taper and HP-only must equal the
    scatter evaluation at the matching cells (exercises assemble_grid)."""
    P = np.array([100.0, 750.0, 3000.0]); T = np.array([280.0, 320.0]); m = np.array([1.0, 3.0])
    G = np.empty(3, object); G[0], G[1], G[2] = P, T, m
    grid = np.asarray(sf.getProp(G, 'NaClaq', DEFPATH, 'rho').rho, float)
    assert grid.shape == (3, 2, 2)
    pts = _scatter(*[(p, t, mm) for p in P for t in T for mm in m])
    scat = np.asarray(sf.getProp(pts, 'NaClaq', DEFPATH, 'rho').rho, float).reshape(3, 2, 2)
    np.testing.assert_allclose(grid, scat, rtol=1e-10)


def test_nacl_zero_molality_runs_and_bulk_props_finite():
    """m = 0 is inside the knot range; bulk properties must be finite there.
    (Partial-molar quantities are ill-conditioned at infinite dilution by
    design and are not asserted.)"""
    out = sf.getProp(_scatter((100.0, 280.0, 0.0)), 'NaClaq', DEFPATH, 'rho', 'G', 'muw', 'aw')
    for k in ('rho', 'G', 'muw', 'aw'):
        assert np.isfinite(getattr(out, k)).all(), k
    np.testing.assert_allclose(out.aw, 1.0, atol=1e-3)
