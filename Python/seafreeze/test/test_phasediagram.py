"""Tests for the full phase diagrams (seafreeze.phasediagram)."""
import numpy as np
import pytest

import seafreeze as sf
from seafreeze import phasediagram as pd

# Literature triple points [T_K, P_MPa] (SF_PhaseLines / phaselines tables)
LIT = {('Ih', 'II', 'III'): (238.237, 209.885), ('II', 'III', 'V'): (249.418, 355.504),
       ('II', 'V', 'VI'): (201.934, 670.840), ('Ih', 'III', 'L'): (251.165, 207.593),
       ('III', 'L', 'V'): (256.164, 350.110), ('L', 'V', 'VI'): (273.407, 634.400)}


@pytest.fixture(scope='module')
def pm():
    return pd.phase_map(np.geomspace(1e-9, 3e3, 500), np.linspace(180, 420, 240))


def test_phase_map_shapes_and_fields(pm):
    assert pm.stable.shape == (500, 240)
    assert pm.names == ['water_Brown2026', 'Ih', 'II', 'III', 'V', 'VI']
    # every phase has a stability field
    for k in range(len(pm.names)):
        assert (pm.stable == k).any(), pm.names[k]


def test_phase_map_spot_checks():
    P = np.array([1e-6, 0.1, 0.1, 300.0, 1000.0])
    T = np.array([300.0, 260.0, 300.0, 230.0, 250.0])
    m = pd.phase_map(P, T)
    names = [m.names[m.stable[i, i]] for i in range(P.size)]
    assert names == ['water_Brown2026', 'Ih', 'water_Brown2026', 'II', 'VI']
    assert m.rho_stable[0, 0] < 1e-3 and m.rho_stable[2, 2] > 990          # vapour, liquid


def test_triple_points_vs_literature(pm):
    tps = pd.triple_points(pm)
    got = {tuple(sorted(tp['labels'])): (tp['T'], tp['P']) for tp in tps}
    for key, (Tl, Pl) in LIT.items():
        k = tuple(sorted(key))
        assert k in got, key
        T, P = got[k]
        assert abs(T - Tl) < 0.1, (key, T, Tl)
        assert abs(P - Pl) < 1.5, (key, P, Pl)
    # Ih - liquid - vapour
    lv = [tp for tp in tps if sorted(tp['labels']) == ['Ih', 'L', 'V']]
    assert len(lv) == 1
    assert abs(lv[0]['T'] - 273.16) < 0.01 and abs(lv[0]['P'] * 1e6 - 611.657) < 0.5


def test_wpd_functions_render():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig = pd.wpd_PT(nP=120, nT=90)
    assert fig.axes and fig.axes[0].get_xscale() == 'log'
    plt.close(fig)
    fig = pd.wpd_rhoT(nP=240, nT=90, nrho=200)
    assert fig.axes
    plt.close(fig)
    fig, ax = plt.subplots()
    pd.wpd_rhoT(ax=ax, rho=(850, 1700), T=(150, 500), P=(1e-10, 1e4), xscale='linear', nP=300, nT=90,
                nrho=200)
    assert ax.get_xscale() == 'linear'
    plt.close(fig)


def test_melt_T_dq2026():
    Tm = pd.melt_T_dq2026(np.array([1e-4, 0.101325, 300.0, 1000.0, 10000.0]))
    assert abs(Tm[0] - 273.16) < 1e-9
    assert abs(Tm[1] - 273.152) < 0.01
    assert 250 < Tm[2] < 260 and 300 < Tm[3] < 310 and 600 < Tm[4] < 800


# --------------------------------------------------------------------------
# the diagrams as data, and the property of the stable phase
# --------------------------------------------------------------------------
def _pts(P, T):
    o = np.empty(len(P), dtype=object)
    for i, pt in enumerate(zip(P, T)):
        o[i] = (float(pt[0]), float(pt[1]))
    return o


@pytest.fixture(scope='module')
def dPT():
    return pd.phase_diagram_PT(P=(1e-8, 3e3), T=(180, 800), nP=240, nT=160)


@pytest.fixture(scope='module')
def drhoT():
    return pd.phase_diagram_rhoT(nP=500, nT=160, nrho=300)


def test_phase_diagram_PT_data(dPT):
    pm = dPT.pm
    assert pm.stable.shape == (240, 160) and dPT.missing.shape == pm.stable.shape
    assert dPT.boundaries and all(len(b[2]) == len(b[3]) for b in dPT.boundaries)
    # saturation: finite on the stable branch, ending at the critical point
    s = dPT.saturation
    ok = np.isfinite(s.P)
    assert ok.sum() > 50 and s.T[ok].max() > 640 and abs(np.nanmax(s.P) - dPT.critical['P']) < 0.5
    # triple points all inside the window
    assert all(1e-8 <= tp['P'] <= 3e3 and 180 <= tp['T'] <= 800 for tp in dPT.triple_points)
    texts = {t[0] for t in dPT.labels}
    assert {'vapour', 'liquid', 'supercritical fluid', 'Ih', 'VI'} <= texts


def test_phase_diagram_PT_matches_wpd(dPT):
    """wpd_PT draws phase_diagram_PT: same stability field."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, pm = pd.wpd_PT(P=(1e-8, 3e3), T=(180, 800), nP=240, nT=160, return_map=True)
    plt.close(fig)
    assert np.array_equal(pm.stable, dPT.pm.stable)


def test_phase_diagram_rhoT_data(drhoT):
    d = drhoT
    assert d.field.shape == (300, 160) and d.two_phase == len(d.pm.names)
    assert (d.field == d.two_phase).any() and (d.field == 0).any()
    # coexisting densities: liquid denser than vapour along the dome
    s = d.saturation
    ok = np.isfinite(s.rho_A) & np.isfinite(s.rho_B)
    assert ok.any() and np.all(s.rho_A[ok] > s.rho_B[ok])
    assert all(len(c[2]) == len(c[3]) == len(c[4]) for c in d.coexistence)
    # no segment spans the vapour-liquid density jump at the triple point
    for i, j, ri, rj, T in d.coexistence:
        for r in (ri, rj):
            with np.errstate(divide='ignore', invalid='ignore'):
                q = np.abs(np.log(r[1:] / r[:-1]))
            assert not np.any(q[np.isfinite(q)] > np.log(20.0))
    # labels re-derived from the field match the ones built with the diagram
    assert sorted(t[0] for t in pd.rhoT_labels(d)) == sorted(t[0] for t in d.labels)
    lin = {t[0]: t[1] for t in pd.rhoT_labels(d, 'linear')}
    assert lin['supercritical fluid'] > 100
    assert any(t[0] == 'L + V' for t in d.labels)


def test_property_map_PT_equals_getProp(dPT):
    m = pd.property_map(dPT, 'rho', 'Cp', 'vel', 'P')
    pm = dPT.pm
    assert np.array_equal(m.phase, pm.stable) and m.coords == 'PT'
    for k, name in enumerate(pm.names):
        ii = np.argwhere((pm.stable == k) & ~m.ideal_gas)
        i, j = ii[len(ii) // 2]
        o = sf.getProp(_pts([pm.P[i]], [pm.T[j]]), name, sf.seafreeze.defpath, 'rho', 'Cp', 'vel')
        for p in ('rho', 'Cp', 'vel'):
            np.testing.assert_allclose(m.values[p][i, j], getattr(o, p)[0], rtol=1e-10, err_msg=name)
        assert m.values['P'][i, j] == pm.P[i]
    assert np.isfinite(m.values['rho'][pm.stable >= 0]).all()


def test_property_map_PT_ideal_gas(dPT):
    m = pd.property_map(dPT, 'rho', 'Cp')
    pm = dPT.pm
    assert m.ideal_gas.any()
    i, j = np.argwhere(m.ideal_gas)[0]
    assert pm.T[j] < 230 and pm.stable[i, j] == 0
    R = 461.5231157
    np.testing.assert_allclose(m.values['rho'][i, j], pm.P[i] * 1e6 / (R * pm.T[j]), rtol=1e-9)


def test_property_map_rhoT(drhoT):
    d = drhoT
    m = pd.property_map(d, 'P', 'Cp', 'rho')
    assert m.coords == 'rhoT' and np.array_equal(m.phase, d.field)
    assert np.isnan(m.values['P'][d.field == d.two_phase]).all()
    # every single-phase point: getProp at the mapped P gives back the density
    for k in range(len(d.pm.names)):
        ii = np.argwhere((d.field == k) & ~m.ideal_gas)
        if not len(ii):
            continue
        i, j = ii[len(ii) // 2]
        o = sf.getProp(_pts([m.values['P'][i, j]], [d.T[j]]), d.pm.names[k], sf.seafreeze.defpath, 'rho', 'Cp')
        np.testing.assert_allclose(o.rho[0], d.rho[i], rtol=1e-3, err_msg=d.pm.names[k])
        np.testing.assert_allclose(o.Cp[0], m.values['Cp'][i, j], rtol=1e-2, err_msg=d.pm.names[k])


def test_property_map_rejects_bad_input(dPT):
    with pytest.raises(ValueError):
        pd.property_map(dPT, 'mus')
    with pytest.raises(TypeError):
        pd.property_map(dPT.pm, 'rho')
    m = pd.property_map(dPT, 'shear')                       # ice-only: NaN in the fluid
    assert np.isnan(m.values['shear'][dPT.pm.stable == 0]).all()
    assert np.isfinite(m.values['shear'][dPT.pm.stable == 1]).all()
