"""1.2 regression tests: spline cache, equal-length grids, wpd with a Helmholtz liquid."""
import os
import shutil
import warnings

import numpy as np
import pytest

import seafreeze as sf
from seafreeze import seafreeze as sfm
from seafreeze.seafreeze import defpath


def _scatter(*cols):
    out = np.empty(len(cols[0]), dtype=object)
    for i, row in enumerate(zip(*cols)):
        out[i] = tuple(float(v) for v in row)
    return out


# --------------------------------------------------------------------------
# spline cache
# --------------------------------------------------------------------------
def test_spline_cache_returns_same_object():
    a = sfm._load_spline(defpath, 'water3')
    b = sfm._load_spline(defpath, 'water3')
    assert a is b                                    # no second file read
    c = sfm._load_spline(defpath, 'Ih')
    assert c is not a and sfm._load_spline(defpath, 'Ih') is c


def test_spline_cache_reloads_modified_file(tmp_path):
    # a private splines/ tree whose Ih file we can touch
    sub = tmp_path / 'ice_Ih'
    sub.mkdir()
    src = os.path.join(defpath, 'ice_Ih', 'ice_Ih.mat')
    dst = sub / 'ice_Ih.mat'
    shutil.copy(src, dst)
    a = sfm._load_spline(str(tmp_path), 'Ih')
    assert sfm._load_spline(str(tmp_path), 'Ih') is a
    st = os.stat(dst)
    os.utime(dst, (st.st_atime, st.st_mtime + 10))  # file changed on disk
    b = sfm._load_spline(str(tmp_path), 'Ih')
    assert b is not a                                # reloaded
    np.testing.assert_array_equal(b['coefs'], a['coefs'])


def test_cached_spline_gives_identical_results():
    pts = _scatter([0.1, 100.0], [300.0, 350.0])
    first = sf.getProp(pts, 'water3', defpath, 'rho', 'Cp')
    second = sf.getProp(pts, 'water3', defpath, 'rho', 'Cp')
    np.testing.assert_array_equal(first.rho, second.rho)
    np.testing.assert_array_equal(first.Cp, second.Cp)


# --------------------------------------------------------------------------
# grids whose axes have the same length (np.array([P, T], dtype=object)
# silently becomes 2-D); used to crash lbftd's Gibbs evaluator
# --------------------------------------------------------------------------
@pytest.mark.parametrize('phase, P, T', [
    ('water1', [100.0, 500.0, 1000.0], [280.0, 300.0, 320.0]),
    ('Ih', [0.1, 50.0, 100.0], [250.0, 260.0, 270.0]),
    ('VI', [800.0, 1000.0, 1200.0], [250.0, 260.0, 270.0]),
])
def test_equal_length_grid_matches_scatter(phase, P, T):
    grid = np.array([np.array(P), np.array(T)], dtype=object)
    assert grid.shape == (2, 3)                      # the problematic 2-D object array
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        g = sf.getProp(grid, phase, defpath, 'G', 'rho', 'Cp')
        Pm, Tm = np.meshgrid(P, T, indexing='ij')
        s = sf.getProp(_scatter(Pm.ravel(), Tm.ravel()), phase, defpath, 'G', 'rho', 'Cp')
    assert g.rho.shape == (3, 3)
    for k in ('G', 'rho', 'Cp'):
        np.testing.assert_allclose(getattr(g, k).ravel(), getattr(s, k), rtol=1e-12)


def test_equal_length_grid_nacl():
    P = np.array([100.0, 200.0, 300.0]); T = np.array([280.0, 290.0, 300.0]); m = np.array([0.5, 1.0, 2.0])
    grid = np.array([P, T, m], dtype=object)
    assert grid.shape == (3, 3)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        g = sf.getProp(grid, 'NaClaq_LP', defpath, 'G', 'rho')
        s = sf.getProp(_scatter([200.0], [290.0], [1.0]), 'NaClaq_LP', defpath, 'G', 'rho')
    assert np.squeeze(g.rho).shape == (3, 3, 3)
    np.testing.assert_allclose(np.squeeze(g.rho)[1, 1, 1], s.rho[0], rtol=1e-12)


# --------------------------------------------------------------------------
# wpd with the Helmholtz liquid
# --------------------------------------------------------------------------
def test_wpd_with_water3_liquid():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        fig = sf.wpd(liquid='water3')
    ax = fig.axes[0]
    melting = [l for l in ax.get_lines() if l.get_label().endswith('-water3')]
    assert len(melting) >= 4                         # Ih, III, V, VI melting curves
    # water3 melting of ice Ih at the lowest plotted pressure is ~273 K
    ih = [l for l in melting if l.get_label() == 'Ih-water3'][0]
    x, y = ih.get_data()
    assert abs(np.asarray(y)[np.nanargmin(np.asarray(x))] - 273.15) < 0.2
    plt.close(fig)


def test_wpd_rejects_unsupported_liquid():
    with pytest.raises(ValueError):
        sf.wpd(liquid='water2')
