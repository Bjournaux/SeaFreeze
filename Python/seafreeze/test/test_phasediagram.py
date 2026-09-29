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
    assert pm.names == ['water3', 'Ih', 'II', 'III', 'V', 'VI']
    # every phase has a stability field
    for k in range(len(pm.names)):
        assert (pm.stable == k).any(), pm.names[k]


def test_phase_map_spot_checks():
    P = np.array([1e-6, 0.1, 0.1, 300.0, 1000.0])
    T = np.array([300.0, 260.0, 300.0, 230.0, 250.0])
    m = pd.phase_map(P, T)
    names = [m.names[m.stable[i, i]] for i in range(P.size)]
    assert names == ['water3', 'Ih', 'water3', 'II', 'VI']
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
