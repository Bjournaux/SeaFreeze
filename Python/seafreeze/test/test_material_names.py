"""1.2 material rename: water1/water2/water3 -> author-year names.

The old names keep working through SeaFreeze 1.x: every public entry point
maps them to the new name with a SeaFreezeDeprecationWarning, shown once per
session per old name, and the results are identical.
"""
import warnings

import numpy as np
import pytest

import seafreeze as sf
from seafreeze import seafreeze as sfm
from seafreeze.seafreeze import defpath

RENAMED = [('water1', 'water_Bollengier2019'), ('water2', 'water_Brown2018'), ('water3', 'water_Brown2026')]
POINT = {'water1': (100.0, 300.0), 'water2': (1000.0, 400.0), 'water3': (0.1, 300.0)}


def _pts(P, T):
    o = np.empty(1, dtype=object)
    o[0] = (float(P), float(T))
    return o


def _record(fn):
    """Run fn with every warning recorded; returns (result, list of deprecation warnings)."""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        out = fn()
    return out, [x for x in w if issubclass(x.category, sf.SeaFreezeDeprecationWarning)]


@pytest.fixture(autouse=True)
def _fresh_session():
    sfm._warned_aliases.clear()          # each test starts as a new session
    yield
    sfm._warned_aliases.clear()


NACL_RENAMED = [('NaClaq_LP', 'NaClaq_Brown2026_LP'), ('NaClaq_HP', 'NaClaq_Brown2026_HP'),
                ('NaClaq_5GPa_2024', 'NaClaq_Brown2024')]
NACL_POINT = {'NaClaq_LP': (100.0, 300.0, 1.0), 'NaClaq_HP': (1500.0, 400.0, 1.0),
              'NaClaq_5GPa_2024': (100.0, 300.0, 1.0)}


def _pts3(P, T, m):
    o = np.empty(1, dtype=object)
    o[0] = (float(P), float(T), float(m))
    return o


def test_alias_table():
    assert sf.MATERIAL_ALIASES == dict(RENAMED + NACL_RENAMED)
    assert sf.MATERIAL_SHORTCUTS == {'NaClaq': 'NaClaq_Brown2026'}
    assert issubclass(sf.SeaFreezeDeprecationWarning, FutureWarning)   # shown by default
    for old, new in RENAMED:
        assert old not in [k for k in sf.phases]                        # listings show new names only
        assert new in [k for k in sf.phases]


@pytest.mark.parametrize('old, new', RENAMED)
def test_getProp_old_name_warns_once_and_matches(old, new):
    pts = _pts(*POINT[old])
    a, w1 = _record(lambda: sf.getProp(pts, old, defpath, 'rho', 'G', 'Cp'))
    b, w2 = _record(lambda: sf.getProp(pts, old, defpath, 'rho', 'G', 'Cp'))
    c, w3 = _record(lambda: sf.getProp(pts, new, defpath, 'rho', 'G', 'Cp'))
    assert len(w1) == 1 and not w2 and not w3                          # once per session, never for new
    msg = str(w1[0].message)
    assert f"'{old}'" in msg and f"'{new}'" in msg and '2.0' in msg
    assert w1[0].filename == __file__                                  # points at the caller's line
    for k in ('rho', 'G', 'Cp'):
        np.testing.assert_array_equal(getattr(a, k), getattr(c, k))
        np.testing.assert_array_equal(getattr(b, k), getattr(c, k))


def test_phases_table_accepts_old_names():
    d, w = _record(lambda: sf.phases['water2'])
    assert d is sf.phases['water_Brown2018'] and len(w) == 1
    got, w = _record(lambda: ('water1' in sf.phases, sf.phases.get('water3')))
    assert got[0] and got[1] is sf.phases['water_Brown2026'] and len(w) == 1
    with pytest.raises(KeyError):
        sf.phases['water4']


def test_old_names_in_the_other_entry_points():
    g = np.empty(2, dtype=object)
    g[0] = np.array([0.1, 300.0, 1000.0]); g[1] = np.array([250.0, 300.0])
    pairs = [
        (lambda: sf.whichphase(g, 'water1'), lambda: sf.whichphase(g, 'water_Bollengier2019')),
        (lambda: sf.phase_lines('Ih', 'water1', segment='stable').T,
         lambda: sf.phase_lines('Ih', 'water_Bollengier2019', segment='stable').T),
        (lambda: sf.phase_range('water2').P, lambda: sf.phase_range('water_Brown2018').P),
        (lambda: sf.rho2P([1000.0], [300.0], 'water1'), lambda: sf.rho2P([1000.0], [300.0], 'water_Bollengier2019')),
        (lambda: sf.saturation([400.0], fluid='water3').P, lambda: sf.saturation([400.0]).P),
        (lambda: sf.sublimation([250.0], fluid='water3').P, lambda: sf.sublimation([250.0]).P),
        (lambda: sf.phase_map([0.1], [300.0], fluid='water3').stable,
         lambda: sf.phase_map([0.1], [300.0]).stable),
    ]
    for old_fn, new_fn in pairs:
        sfm._warned_aliases.clear()
        a, w = _record(old_fn)
        b, wn = _record(new_fn)
        assert len(w) == 1 and not wn
        np.testing.assert_array_equal(np.asarray(a, float), np.asarray(b, float))


def test_phasenum2phase_returns_new_names():
    assert sf.phasenum2phase(0) == 'water_Bollengier2019'
    out, w = _record(lambda: sf.phasenum2phase(0, 'water3'))
    assert out == 'water_Brown2026' and len(w) == 1
    assert sf.phasenum2phase(6) == 'VI'


def test_unknown_names_still_fail():
    with pytest.raises(ValueError):
        sf.getProp(_pts(0.1, 300.0), 'water4', defpath)


@pytest.mark.parametrize('old, new', NACL_RENAMED)
def test_nacl_old_names_warn_once_and_match(old, new):
    pts = _pts3(*NACL_POINT[old])
    a, w1 = _record(lambda: sf.getProp(pts, old, defpath, 'rho', 'muw'))
    b, w2 = _record(lambda: sf.getProp(pts, old, defpath, 'rho', 'muw'))
    c, w3 = _record(lambda: sf.getProp(pts, new, defpath, 'rho', 'muw'))
    assert len(w1) == 1 and not w2 and not w3
    assert f"'{old}'" in str(w1[0].message) and f"'{new}'" in str(w1[0].message)
    for k in ('rho', 'muw'):
        np.testing.assert_array_equal(getattr(a, k), getattr(c, k))


def test_naclaq_shortcut_is_silent_and_identical():
    pts = _pts3(200.0, 290.0, 1.5)
    a, w = _record(lambda: sf.getProp(pts, 'NaClaq', defpath, 'rho', 'muw'))
    b, _ = _record(lambda: sf.getProp(pts, 'NaClaq_Brown2026', defpath, 'rho', 'muw'))
    assert not w                                          # a shortcut, not a deprecation
    np.testing.assert_array_equal(a.rho, b.rho)
    np.testing.assert_array_equal(a.muw, b.muw)
    assert 'NaClaq' in sf.phases and sf.phases['NaClaq'] is sf.phases['NaClaq_Brown2026']
    r, w = _record(lambda: sf.phase_lines('Ih', 'NaClaq', m=1.0, segment='stable'))
    assert not w and r.matB == 'NaClaq_Brown2026'
    g = np.empty(3, dtype=object)
    g[0] = np.array([0.1, 100.0]); g[1] = np.array([260.0, 300.0]); g[2] = np.array([1.0])
    np.testing.assert_array_equal(sf.whichphase(g, 'NaClaq'), sf.whichphase(g, 'NaClaq_Brown2026'))
