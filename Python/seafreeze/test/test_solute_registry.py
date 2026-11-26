"""Solute registry: chemistry functions key on the registry, not on the name 'NaClaq'.

The central check registers a *second* aqueous binary at test time, backed by
existing NaCl spline files under new material names, and requires every
chemistry-aware function to treat it exactly like the built-in material it
aliases.  If that holds, adding a real binary needs a registry entry and a
spline file, and no change to getProp / whichphase / phase_lines / wpd / rho2P.
"""
import pathlib
import re
import warnings

import numpy as np
import pytest

import seafreeze as sf
import seafreeze.seafreeze as sfm
from seafreeze import phaselines as pl

MW_NACL = sfm.mNaCl_kgmol
SRC = pathlib.Path(sfm.__file__).parent


def _scatter(pts):
    PTm = np.empty((len(pts),), dtype=object)
    for i, p in enumerate(pts):
        PTm[i] = tuple(float(v) for v in p)
    return PTm


def _grid(P, T, m):
    return np.array([np.asarray(P, float), np.asarray(T, float), np.asarray(m, float)], dtype=object)


def _assert_same_props(a, b, label):
    assert set(vars(a)) == set(vars(b)), f'{label}: property sets differ'
    for k, v in vars(a).items():
        np.testing.assert_array_equal(np.asarray(v, float), np.asarray(getattr(b, k), float),
                                      err_msg=f'{label}: {k}')


@pytest.fixture
def dummy():
    """'TESTaq' aliases NaClaq_5GPa_2024 (single spline); 'TESTst' aliases stitched NaClaq."""
    sf.register_solution('TESTaq', 'TEST', MW_NACL, 2, 0.0002,
                         spline=sfm._SPLINE_MAP['NaClaq_5GPa_2024'])
    sf.register_solution('TESTlp', 'TEST', MW_NACL, 2, 0.0002, spline=sfm._SPLINE_MAP['NaClaq_LP'])
    sf.register_solution('TESThp', 'TEST', MW_NACL, 2, 0.0002, spline=sfm._SPLINE_MAP['NaClaq_HP'])
    sf.register_solution('TESTst', 'TEST', MW_NACL, 2, 0.0002,
                         stitch=('TESTlp', 'TESThp', sfm._STITCH_P_LO, sfm._STITCH_P_HI))
    yield
    for name in ('TESTst', 'TESThp', 'TESTlp', 'TESTaq'):
        sf.unregister_solution(name)


# ---------------------------------------------------------------------------
# Registry queries
# ---------------------------------------------------------------------------
def test_builtin_classification():
    for mat in ('NaClaq', 'NaClaq_LP', 'NaClaq_HP', 'NaClaq_5GPa_2024'):
        assert sf.is_solution(mat), mat
        assert sf.phases[mat].solute == 'NaCl'
        assert not sf.phases[mat].experimental
    assert sf.is_solution('NH3aq_RK') and sf.phases['NH3aq_RK'].solute == 'NH3' and sf.phases['NH3aq_RK'].nu == 1
    assert sf.phases['NH3aq_RK'].experimental
    for mat in ('Ih', 'II', 'III', 'V', 'VI', 'VII_X_French', 'water1', 'water2', 'water_IAPWS95'):
        assert not sf.is_solution(mat), mat
    assert not sf.is_solution('NotAPhase')
    assert sfm.is_stitched('NaClaq') and not sfm.is_stitched('NaClaq_LP')


def test_solute_lookup():
    assert sf.solution_phases('NaCl') == ['NaClaq', 'NaClaq_LP', 'NaClaq_HP', 'NaClaq_5GPa_2024']
    assert sf.solution_phases('nacl') == sf.solution_phases('NaCl')
    assert sf.solution_phases('NH3') == ['NH3aq_RK'] and sf.solute_phase('nh3') == 'NH3aq_RK'
    assert sf.solution_phases() == sf.solution_phases('NaCl') + sf.solution_phases('NH3')
    assert sf.solution_phases('MgSO4') == []
    assert sf.solute_phase('NaCl') == 'NaClaq'
    assert sf.solute_phase('nacl') == 'NaClaq'
    assert sf.solute_phase('NaClaq_LP') == 'NaClaq_LP'
    with pytest.raises(ValueError, match='Unknown solute'):
        sf.solute_phase('brine')


def test_register_validation(dummy):
    with pytest.raises(ValueError, match='already registered'):
        sf.register_solution('TESTaq', 'TEST', MW_NACL, 2, 0.0002, spline=('a', 'b.mat'))
    with pytest.raises(ValueError, match='already registered'):
        sf.register_solution('Ih', 'TEST', MW_NACL, 2, 0.0002, spline=('a', 'b.mat'))
    with pytest.raises(ValueError, match='exactly one'):
        sf.register_solution('TESTx', 'TEST', MW_NACL, 2, 0.0002)
    with pytest.raises(ValueError, match='exactly one'):
        sf.register_solution('TESTx', 'TEST', MW_NACL, 2, 0.0002, spline=('a', 'b.mat'),
                             stitch=('TESTlp', 'TESThp', 1.0, 2.0))
    with pytest.raises(ValueError, match='Stitch component'):
        sf.register_solution('TESTx', 'TEST', MW_NACL, 2, 0.0002, stitch=('water1', 'TESThp', 1.0, 2.0))
    with pytest.raises(ValueError, match='not a registered solution'):
        sf.unregister_solution('water1')
    assert 'TESTx' not in sf.phases


def test_unregister_restores_registry():
    before = dict(sf.phases)
    sf.register_solution('TESTtmp', 'TEST', MW_NACL, 2, 0.0002, spline=sfm._SPLINE_MAP['NaClaq_LP'])
    assert sf.is_solution('TESTtmp') and sf.solute_phase('TEST') == 'TESTtmp'
    sf.unregister_solution('TESTtmp')
    assert dict(sf.phases) == before and 'TESTtmp' not in sfm._SPLINE_MAP


# ---------------------------------------------------------------------------
# A second binary behaves exactly like the material it aliases
# ---------------------------------------------------------------------------
SCATTER = _scatter([(0.1, 300.0, 0.5), (200.0, 280.0, 2.0), (450.0, 320.0, 4.0), (700.0, 300.0, 1.0),
                    (1500.0, 350.0, 3.0)])
GRID = _grid([0.1, 300.0, 700.0, 1200.0], [270.0, 300.0, 340.0], [0.0, 1.0, 3.0])


@pytest.mark.parametrize('alias, builtin', [('TESTaq', 'NaClaq_5GPa_2024'), ('TESTst', 'NaClaq')])
def test_getProp_identical(dummy, alias, builtin):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        _assert_same_props(sf.getProp(SCATTER, alias), sf.getProp(SCATTER, builtin), f'{alias} scatter')
        _assert_same_props(sf.getProp(GRID, alias), sf.getProp(GRID, builtin), f'{alias} grid')
        sel = ('rho', 'muw', 'Vw', 'xs', 'm', 'Js')
        _assert_same_props(sf.getProp(SCATTER, alias, sfm.defpath, *sel),
                           sf.getProp(SCATTER, builtin, sfm.defpath, *sel), f'{alias} selective')


@pytest.mark.parametrize('alias, builtin', [('TESTaq', 'NaClaq_5GPa_2024'), ('TESTst', 'NaClaq')])
def test_whichphase_identical(dummy, alias, builtin):
    pts = _scatter([(100.0, 250.0, 3.0), (900.0, 285.0, 3.0), (300.0, 240.0, 1.0), (650.0, 260.0, 2.0)])
    np.testing.assert_array_equal(sf.whichphase(pts, alias), sf.whichphase(pts, builtin))
    grid = _grid(np.arange(0.1, 900.0, 100.0), np.arange(235.0, 300.0, 10.0), [1.0])
    np.testing.assert_array_equal(sf.whichphase(grid, alias), sf.whichphase(grid, builtin))


def test_whichphase_accepts_solute_name(dummy):
    pts = _scatter([(100.0, 250.0, 3.0), (900.0, 285.0, 3.0)])
    np.testing.assert_array_equal(sf.whichphase(pts, 'NaCl'), sf.whichphase(pts, 'NaClaq'))
    np.testing.assert_array_equal(sf.whichphase(pts, 'TEST'), sf.whichphase(pts, 'TESTaq'))
    with pytest.raises(ValueError):
        sf.whichphase(pts, 'brine')


def test_phase_range_identical(dummy):
    assert sf.phase_range('TESTaq') == sf.phase_range('NaClaq_5GPa_2024')
    assert sf.phase_range('TESTst') == sf.phase_range('NaClaq')


@pytest.mark.parametrize('ice', ['Ih', 'III', 'V', 'VI', 'II'])
def test_phase_lines_identical_to_NaClaq(dummy, ice):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        a = sf.phase_lines(ice, 'TESTst', m=1.0)
        b = sf.phase_lines(ice, 'NaClaq', m=1.0)
    assert (a.matA, a.matB) == (ice, 'TESTst') and a.m == 1.0
    np.testing.assert_array_equal(a.P, b.P)
    np.testing.assert_array_equal(a.T, b.T)
    np.testing.assert_array_equal(a.stable, b.stable)
    np.testing.assert_array_equal(a.triple_points, b.triple_points)


def test_phase_lines_argument_order_and_vector_m(dummy):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        a = sf.phase_lines('TESTst', 'Ih', m=1.0)
        b = sf.phase_lines('Ih', 'TESTst', m=1.0)
        both = sf.phase_lines('Ih', 'TESTst', m=[0.5, 2.0])
    np.testing.assert_array_equal(a.T, b.T)
    assert [r.m for r in both] == [0.5, 2.0]


def test_phase_lines_validation_names_the_solution(dummy):
    with pytest.raises(ValueError, match='involves TESTst'):
        sf.phase_lines('Ih', 'TESTst')
    with pytest.raises(ValueError, match='outside TESTst spline range'):
        sf.phase_lines('Ih', 'TESTst', m=100.0)
    with pytest.raises(ValueError):
        sf.phase_lines('water1', 'TESTst', m=1.0)       # liquid-liquid pair is not defined
    with pytest.raises(ValueError):
        sf.phase_lines('Ih', 'water1', m=0.5)


def test_single_spline_variant_gets_melting_lines():
    """Any registered solution has ice melting lines, not only the default 'NaClaq'."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        lp = sf.phase_lines('Ih', 'NaClaq_LP', m=1.0)
        st = sf.phase_lines('Ih', 'NaClaq', m=1.0)
    # Ih melts below 210 MPa, where stitched NaClaq is the LP spline alone
    np.testing.assert_allclose(np.interp([50.0, 150.0], lp.P[np.argsort(lp.P)], lp.T[np.argsort(lp.P)]),
                               np.interp([50.0, 150.0], st.P[np.argsort(st.P)], st.T[np.argsort(st.P)]),
                               atol=0.05)


def test_wpd_overlay_by_solute_or_material(dummy):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        base = len(sf.wpd().axes[0].get_lines())
        by_solute = sf.wpd(solute='TEST', m=[1.0])           # resolves to first registered: TESTaq
        by_material = sf.wpd(solute='TESTst', m=[1.0])
        labels = [l.get_label() for l in by_material.axes[0].get_lines()]
        assert len(by_solute.axes[0].get_lines()) > base
        assert any(lab.startswith('Ih-TESTst') for lab in labels)
        with pytest.raises(ValueError, match='Unknown solute'):
            sf.wpd(solute='brine', m=[1.0])
        with pytest.raises(ValueError):
            sf.wpd(solute='TEST')
    plt.close('all')


def test_rho2P_round_trip(dummy):
    P_ref = np.array([50.0, 400.0, 800.0])
    PTm = _scatter([(p, 300.0, 2.0) for p in P_ref])
    for alias in ('TESTaq', 'TESTst'):
        rho = sf.getProp(PTm, alias, sfm.defpath, 'rho').rho
        np.testing.assert_allclose(sf.rho2P(rho, 300.0, alias, m=2.0), P_ref, atol=0.05, err_msg=alias)
        with pytest.raises(ValueError, match='requires the'):
            sf.rho2P(rho, 300.0, alias)


def test_load_stitched_material_names_its_components(dummy):
    with pytest.raises(ValueError, match="'TESTlp' / 'TESThp'"):
        sfm._load_spline(sfm.defpath, 'TESTst')


# ---------------------------------------------------------------------------
# Guard: the name is not hard-coded in the chemistry code paths any more
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('fname', ['seafreeze.py', 'phaselines.py', 'rho2P.py'])
def test_no_name_based_dispatch(fname):
    code = (SRC / fname).read_text()
    # strip comments and docstrings
    code = re.sub(r'"""(?:.|\n)*?"""', '', code)
    code = re.sub(r'#.*', '', code)
    assert "startswith('NaClaq')" not in code and 'startswith("NaClaq")' not in code
    for m in re.finditer(r"(==|!=|\bin\b)\s*\(?['\"]NaClaq", code):
        ctx = code[max(0, m.start() - 60):m.end() + 20]
        assert '_SOLUTION_PAIR_TEMPLATE' in ctx, f'{fname}: name-based dispatch remains: {ctx!r}'
