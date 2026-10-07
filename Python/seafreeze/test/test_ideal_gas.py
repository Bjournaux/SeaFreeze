"""lbftd.evalHelmholtz.ideal_gas_props: the psi surface's ideal-gas part alone.

The dilute-vapour extension below the water_Brown2026 surface (sublimation,
phase_map); checked against the full surface where both exist, against
thermodynamic identities and finite differences everywhere else.
"""
import warnings

import numpy as np
import pytest

import seafreeze as sf
from seafreeze.seafreeze import defpath, _load_spline
from lbftd import evalHelmholtz as eh


@pytest.fixture(scope='module')
def sp3():
    return _load_spline(defpath, 'water_Brown2026')


@pytest.fixture(scope='module')
def R(sp3):
    return eh._scalar(sp3, 'R')


def _scatter(a, b):
    out = np.empty(len(a), dtype=object)
    for i, (x, y) in enumerate(zip(a, b)):
        out[i] = (float(x), float(y))
    return out


def _relerr(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    return np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-300))


def test_same_fields_as_getProp_water3(sp3):
    o = eh.ideal_gas_props(sp3, 1e-7, 300.0)
    w = sf.getProp(_scatter([0.1], [300.0]), 'water_Brown2026')
    assert set(vars(o)) == set(vars(w))


def test_matches_water3_in_the_dilute_limit(sp3):
    # the real fluid tends to its ideal-gas part as P -> 0.  Above ~700 K
    # the surface's reacting-mixture reference term takes over, so stay below.
    T = np.repeat(np.arange(300.0, 651.0, 50.0), 3)
    P = np.tile([1e-7, 1e-8, 1e-9], T.size // 3)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        w = sf.getProp(_scatter(P, T), 'water_Brown2026')
    o = eh.ideal_gas_props(sp3, P, T)
    for k in ('G', 'S', 'H', 'U', 'A', 'rho', 'Cp', 'Cv', 'vel', 'Kt', 'Ks', 'alpha',
              'Js', 'gamma_Gruneisen'):
        assert _relerr(getattr(o, k), getattr(w, k)) < 1e-5, k
    np.testing.assert_array_equal(o.P, P)
    np.testing.assert_array_equal(o.T, T)


def test_identities(sp3, R):
    T = np.array([100.0, 150.0, 229.0, 273.16, 500.0, 1000.0, 2000.0])
    P = np.array([1e-12, 1e-9, 1e-6, 6e-4, 1e-3, 0.1, 1.0])
    o = eh.ideal_gas_props(sp3, P, T)
    assert _relerr(o.Cp - o.Cv, np.full(T.size, R)) < 1e-12
    assert _relerr(o.Ks / o.Kt, o.Cp / o.Cv) < 1e-12
    assert _relerr(o.vel ** 2, o.Ks * 1e6 / o.rho) < 1e-12
    assert _relerr(o.P * 1e6, o.rho * R * T) < 1e-12              # Z = 1
    assert _relerr(o.Kt, P) < 1e-15 and np.all(o.Kp == 1.0)
    assert _relerr(o.alpha * T, np.ones(T.size)) < 1e-15
    assert _relerr(o.G, o.H - T * o.S) < 1e-10
    assert _relerr(o.A, o.U - T * o.S) < 1e-10
    assert _relerr(o.G - o.A, R * T) < 1e-10
    # derived props with getProp's definitions (Js: isentropic dT/dP, K/MPa)
    assert _relerr(o.Js, T * o.alpha / (o.rho * o.Cp) * 1e6) < 1e-12
    assert _relerr(o.gamma_Gruneisen, o.alpha * o.Kt * 1e6 / (o.rho * o.Cv)) < 1e-12


def test_finite_differences(sp3):
    # S = -(dA/dT)_rho and Cv = T (dS/dT)_rho, including below the surface
    rho = np.array([1e-9, 1e-6, 1e-3, 1e-1])
    T = np.array([150.0, 200.0, 400.0, 900.0])
    h = 1e-4
    o = eh.ideal_gas_props(sp3, rho, T, rhoT=True)
    up = eh.ideal_gas_props(sp3, rho, T * (1 + h), 'A', 'S', rhoT=True)
    dn = eh.ideal_gas_props(sp3, rho, T * (1 - h), 'A', 'S', rhoT=True)
    assert _relerr(-(up.A - dn.A) / (2 * h * T), o.S) < 1e-7
    assert _relerr(T * (up.S - dn.S) / (2 * h * T), o.Cv) < 1e-7


def test_rhoT_and_P_input_agree(sp3):
    P = np.array([1e-10, 1e-7, 1e-4, 1e-2]); T = np.array([150.0, 250.0, 300.0, 800.0])
    a = eh.ideal_gas_props(sp3, P, T)
    b = eh.ideal_gas_props(sp3, a.rho, T, rhoT=True)
    for k, v in vars(a).items():
        np.testing.assert_allclose(getattr(b, k), v, rtol=1e-13, err_msg=k)


def test_finite_below_the_surface(sp3):
    Tmin = eh.domain(sp3)[1][0]
    assert Tmin > 150                                    # the surface stops at 230 K
    P = np.geomspace(1e-14, 1e-6, 9)
    o = eh.ideal_gas_props(sp3, P, 150.0)
    for k, v in vars(o).items():
        assert v.shape == P.shape and np.all(np.isfinite(v)), k
    assert np.all(o.Cv > 0) and np.all(o.S > 0) and np.all(o.vel > 0)


def test_G_rho_equal_ideal_gas(sp3):
    P = np.array([1e-9, 1e-7, 6e-4]); T = np.array([150.0, 220.0, 273.16])
    G, rho = eh.ideal_gas(sp3, P, T)
    o = eh.ideal_gas_props(sp3, P, T, 'G', 'rho')
    np.testing.assert_array_equal(o.G, G)
    np.testing.assert_array_equal(o.rho, rho)


def test_broadcasting_and_selection(sp3):
    P = np.geomspace(1e-9, 1e-5, 4)[:, None]; T = np.array([150.0, 200.0, 250.0])
    o = eh.ideal_gas_props(sp3, P, T, 'Cp', 'P', 'T')
    assert set(vars(o)) == {'Cp', 'P', 'T'}
    assert o.Cp.shape == o.P.shape == o.T.shape == (4, 3)
    s = eh.ideal_gas_props(sp3, 1e-7, 200.0)
    assert np.ndim(s.G) == 0 and np.isfinite(s.G)
    with pytest.raises(ValueError):
        eh.ideal_gas_props(sp3, 1e-7, 200.0, 'shear')
    with pytest.raises(ValueError):
        eh.ideal_gas_props({'eos': 'F_rhoT'}, 1e-7, 200.0)
