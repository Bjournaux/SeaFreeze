"""Grid output vs scatter evaluation, point by point, for (P,T) and (rho,T) grids.

Regression for getProp(..., rhoT=True) on grids: the density echo is the 1-D
input axis, and the derived Js / gamma_Gruneisen once used it as if it were
the full grid -- a broadcast error on non-square grids and silently wrong
values (rho varying along the T axis) on square ones.
"""
import warnings

import numpy as np
import pytest

import seafreeze as sf
from seafreeze.seafreeze import defpath


def _grid(*axes):
    g = np.empty(len(axes), dtype=object)
    for i, a in enumerate(axes):
        g[i] = np.asarray(a, float)
    return g


def _scatter(*cols):
    out = np.empty(len(cols[0]), dtype=object)
    for i, row in enumerate(zip(*cols)):
        out[i] = tuple(float(v) for v in row)
    return out


def _get(PTm, phase, *props, **kw):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return sf.getProp(PTm, phase, defpath, *props, **kw)


def _assert_grid_matches_scatter(g, s, axes, rhoT):
    """Every computed field of the grid result has the full grid shape and
    equals the scatter result at the same points.  Coordinate echoes are the
    input axes: 1-D on a 2-D grid, broadcastable (n,1,1)-style on a 3-D one."""
    shape = tuple(a.size for a in axes)
    mesh = np.meshgrid(*axes, indexing='ij')
    echo = ('rho', 'T') if rhoT else ('P', 'T', 'm')
    for k, v in vars(g).items():
        v = np.asarray(v, float)
        if k in echo:
            i = echo.index(k)
            if len(shape) == 2:
                np.testing.assert_array_equal(v, axes[i], err_msg=k)
            else:
                np.testing.assert_array_equal(np.broadcast_to(v, shape), mesh[i], err_msg=k)
            continue
        assert v.shape == shape, (k, v.shape)
        ref = np.asarray(getattr(s, k), float)
        np.testing.assert_array_equal(np.isnan(v.ravel()), np.isnan(ref), err_msg=k)
        np.testing.assert_allclose(v.ravel(), ref, rtol=1e-9, atol=0, equal_nan=True, err_msg=k)


# (rho, T) axes: non-square and square.  water_Brown2026 spans the dilute vapour to
# the dense liquid; the square case exercises the silent broadcast.
_W3_NONSQ = (np.r_[np.geomspace(1e-5, 10, 5), np.linspace(900, 1300, 6)], np.linspace(250, 700, 4))
_W3_SQ = (np.linspace(950, 1150, 4), np.array([260., 300, 400, 500]))


@pytest.mark.parametrize('axes', [_W3_NONSQ, _W3_SQ], ids=['nonsquare', 'square'])
def test_water3_rhoT_grid_matches_scatter(axes):
    g = _get(_grid(*axes), 'water_Brown2026', rhoT=True)
    R, T = np.meshgrid(*axes, indexing='ij')
    s = _get(_scatter(R.ravel(), T.ravel()), 'water_Brown2026', rhoT=True)
    assert g.P.shape == R.shape and g.Js.shape == R.shape
    assert np.array_equal(g.rho, axes[0]) and np.array_equal(g.T, axes[1])
    _assert_grid_matches_scatter(g, s, axes, rhoT=True)
    # Js and gamma really are per point (the square grid used to mix axes)
    assert np.isfinite(g.Js).sum() > 0.5 * g.Js.size


def test_water3_rhoT_large_nonsquare_grid():
    # the case that crashed: (800, 420) grid, vapour to compressed liquid
    rho = np.geomspace(1e-7, 4e3, 80); T = np.linspace(150, 1800, 42)
    g = _get(_grid(rho, T), 'water_Brown2026', rhoT=True)
    for k, v in vars(g).items():
        assert np.shape(v) == ((80,) if k == 'rho' else (42,) if k == 'T' else (80, 42)), k
    i = np.array([5, 40, 60, 75]); j = np.array([10, 3, 20, 41])
    s = _get(_scatter(rho[i], T[j]), 'water_Brown2026', 'Js', 'gamma_Gruneisen', 'P', rhoT=True)
    for k in ('Js', 'gamma_Gruneisen', 'P'):
        np.testing.assert_allclose(getattr(g, k)[i, j], getattr(s, k), rtol=1e-9, equal_nan=True, err_msg=k)


@pytest.mark.parametrize('props', [('Js',), ('gamma_Gruneisen',), ('rho',), ('P',), ('T',),
                                   ('Js', 'P'), ('rho', 'Js', 'T')])
@pytest.mark.parametrize('phase', ['water_Brown2026', 'water_Bollengier2019'])
def test_rhoT_grid_selective_props(phase, props):
    axes = (np.linspace(1000, 1100, 5), np.linspace(280, 330, 3)) if phase == 'water_Bollengier2019' else _W3_NONSQ
    g = _get(_grid(*axes), phase, *props, rhoT=True)
    assert set(vars(g)) == set(props)                # nothing leaks (no bare prerequisites)
    R, T = np.meshgrid(*axes, indexing='ij')
    s = _get(_scatter(R.ravel(), T.ravel()), phase, *props, rhoT=True)
    _assert_grid_matches_scatter(g, s, axes, rhoT=True)


@pytest.mark.parametrize('axes', [(np.linspace(1000, 1100, 5), np.linspace(280, 330, 3)),
                                  (np.linspace(1000, 1100, 3), np.linspace(280, 330, 3))],
                         ids=['nonsquare', 'square'])
def test_water1_rhoT_grid_matches_scatter(axes):
    g = _get(_grid(*axes), 'water_Bollengier2019', rhoT=True)
    R, T = np.meshgrid(*axes, indexing='ij')
    s = _get(_scatter(R.ravel(), T.ravel()), 'water_Bollengier2019', rhoT=True)
    _assert_grid_matches_scatter(g, s, axes, rhoT=True)


@pytest.mark.parametrize('phase', ['water_Brown2026', 'water_Bollengier2019', 'VI'])
def test_PT_grid_matches_scatter(phase):
    axes = {'water_Brown2026': (np.array([1e-4, 0.1, 50, 500, 1500]), np.array([260., 300, 450])),
            'water_Bollengier2019': (np.array([0.1, 50, 500, 1500]), np.array([260., 300, 350])),
            'VI': (np.array([800., 1000, 1200, 1500]), np.array([250., 270]))}[phase]
    g = _get(_grid(*axes), phase)
    P, T = np.meshgrid(*axes, indexing='ij')
    s = _get(_scatter(P.ravel(), T.ravel()), phase)
    _assert_grid_matches_scatter(g, s, axes, rhoT=False)


@pytest.mark.parametrize('props', [('P',), ('T',), ('Js',), ('Vp',), ('rho', 'gamma_Gruneisen')])
def test_PT_selective_props_strip_prerequisites(props):
    g = _get(_grid([800., 1000, 1200], [250., 270]), 'VI', *props)
    out = set(vars(g))
    assert set(props) <= out
    # prerequisites getProp adds internally never leak (lbftd's own
    # by-products, e.g. vel next to Ks, may)
    assert not (out - set(props)) & {'rho', 'alpha', 'Cp', 'Cv', 'Kt', 'Ks', 'muw', 'G', 'S'}
    # echo-only requests no longer evaluate every property
    if props in (('P',), ('T',)):
        assert out == set(props)


@pytest.mark.parametrize('axes', [(np.array([1000., 1050, 1100, 1150]), np.array([280., 320]), np.array([0.5, 1.0, 2.0])),
                                  (np.array([1000., 1050, 1100]), np.array([280., 300, 320]), np.array([0.5, 1.0, 2.0]))],
                         ids=['nonsquare', 'cube'])
def test_nacl_rhoT_grid_matches_scatter(axes):
    props = ('P', 'rho', 'T', 'G', 'Js', 'gamma_Gruneisen', 'muw')
    g = _get(_grid(*axes), 'NaClaq', *props, rhoT=True)
    R, T, M = np.meshgrid(*axes, indexing='ij')
    s = _get(_scatter(R.ravel(), T.ravel(), M.ravel()), 'NaClaq', *props, rhoT=True)
    # echoes broadcast like P and T on a (P,T,m) grid
    assert g.rho.shape == (axes[0].size, 1, 1) and g.T.shape == (1, axes[1].size, 1)
    _assert_grid_matches_scatter(g, s, axes, rhoT=True)


@pytest.mark.parametrize('props', [('P',), ('xs',), ('m', 'f'), ('Js',)])
def test_nacl_PT_selective_props(props):
    axes = (np.array([10., 100, 300, 450]), np.array([280., 320]), np.array([0.5, 1.0, 2.0]))
    for phase in ('NaClaq_LP', 'NaClaq'):
        g = _get(_grid(*axes), phase, *props)
        assert set(vars(g)) == set(props), phase
        for k in props:
            assert np.broadcast_to(getattr(g, k), (4, 2, 3)).shape == (4, 2, 3)
