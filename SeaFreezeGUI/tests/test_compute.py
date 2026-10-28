"""Tests for SeaFreezeGUI/core/compute.py.

Calls the uncached function bodies (``.__wrapped__``) so no Streamlit
server is needed.  Inputs mirror exactly what app.py sends: tuples of
floats, a tuple of property symbols, and 'scatter' or 'grid'.
"""

import logging

import numpy as np
import pytest

from core.compute import (compute_properties, get_phase_range,
                          get_stability_boundaries, get_phase_line,
                          get_phase_line_full, _canonical)
from core.constants import (ALL_MATERIALS, available_properties,
                            categorized_properties, is_nacl, is_solid)
import seafreeze as sf
from seafreeze.seafreeze import defpath as DEFPATH

logging.getLogger('streamlit').setLevel(logging.ERROR)

_compute = compute_properties.__wrapped__
_range = get_phase_range.__wrapped__
_bounds = get_stability_boundaries.__wrapped__
_line = get_phase_line.__wrapped__
_line_full = get_phase_line_full.__wrapped__

P1, T1 = (300.0,), (260.0,)
PG, TG = (100.0, 300.0, 500.0), (250.0, 260.0, 270.0)


def _m(material, grid=False):
    if not is_nacl(material):
        return None
    return (0.5, 1.0) if grid else (1.0,)


def test_gui_material_list_matches_package():
    assert set(ALL_MATERIALS) == set(sf.phases)


@pytest.mark.parametrize('material', ALL_MATERIALS)
def test_gui_property_symbols_exist_in_getprop(material):
    """Every symbol the GUI offers must be a property getProp can return."""
    props = tuple(available_properties(material))
    PTm = np.empty(1, object)
    PTm[0] = (1000.0, 300.0, 1.0) if is_nacl(material) else (1000.0, 300.0)
    # pick an in-domain point per material class
    if material == 'VII_X_French':
        PTm[0] = (5000.0, 300.0)
    elif material == 'NaClaq_HP':
        PTm[0] = (2000.0, 600.0, 1.0)
    out = sf.getProp(PTm, material, DEFPATH, *props)
    assert set(vars(out)) == set(props)


@pytest.mark.parametrize('material', ALL_MATERIALS)
def test_categorized_properties_cover_available(material):
    flat = {s for _, items in categorized_properties(material) for s, _, _ in items}
    assert flat == set(available_properties(material))


@pytest.mark.parametrize('material', ALL_MATERIALS)
@pytest.mark.parametrize('mode', ['scatter', 'grid'])
def test_compute_properties_all_selected(material, mode):
    props = tuple(available_properties(material))
    P, T = (P1, T1) if mode == 'scatter' else (PG, TG)
    res = _compute(P, T, _m(material, mode == 'grid'), material, props, mode)
    assert set(res) == set(props), f'missing={set(props) - set(res)} extra={set(res) - set(props)}'
    exp = (1,) if mode == 'scatter' else ((3, 3, 2) if is_nacl(material) else (3, 3))
    for k, v in res.items():
        assert isinstance(v, np.ndarray), k
        if k not in ('xs', 'xw', 'f'):
            assert np.squeeze(v).shape in (exp, ()), f'{k}: {v.shape}'


@pytest.mark.parametrize('material', ALL_MATERIALS)
def test_compute_properties_single_selected(material):
    for prop in available_properties(material):
        res = _compute(P1, T1, _m(material), material, (prop,), 'scatter')
        assert set(res) == {prop}


@pytest.mark.parametrize('material', ['VI', 'water1', 'NaClaq'])
def test_compute_properties_selective_matches_all(material):
    sel = _compute(PG, TG, _m(material, True), material, ('rho', 'Cp', 'Js'), 'grid')
    full = _compute(PG, TG, _m(material, True), material, (), 'grid')
    for k in sel:
        np.testing.assert_allclose(sel[k], full[k], rtol=1e-12, equal_nan=True, err_msg=k)


def test_compute_properties_empty_props_returns_everything():
    res = _compute(P1, T1, None, 'VI', (), 'scatter')
    assert set(available_properties('VI')) <= set(res)


def test_compute_properties_out_of_domain_point_is_nan_not_error():
    res = _compute((300.0,), (260.0,), None, 'VII_X_French', ('rho', 'G'), 'scatter')
    assert set(res) == {'rho', 'G'}
    assert np.isnan(res['rho']).all() and np.isnan(res['G']).all()


def test_compute_properties_grid_equal_length_axes():
    """Equal-length P and T axes used to collapse into a 2-D object array."""
    res = _compute((100.0, 200.0, 300.0), (250.0, 260.0, 270.0), None, 'Ih', ('rho',), 'grid')
    assert res['rho'].shape == (3, 3)


@pytest.mark.parametrize('material', ALL_MATERIALS)
def test_get_phase_range(material):
    P, T, m = _range(material)
    assert P[0] < P[1] and T[0] < T[1]
    assert (m is None) == (not is_nacl(material))


def test_canonical_mapping():
    assert _canonical('water2') == 'water1'
    assert _canonical('NaClaq_HP') == 'NaClaq'
    assert _canonical('VI') == 'VI'


def test_stability_boundaries_for_ice_VI():
    b = _bounds('VI')
    pairs = {(a, o) for a, o, _, _ in b}
    assert ('VI', 'water1') in pairs and ('VI', 'V') in pairs
    for _, _, P, T in b:
        assert len(P) == len(T) > 0


def test_phase_line_helpers():
    assert _line('Ih', 'water1') is not None
    full = _line_full('Ih', 'water1')
    assert full is not None
    P, T, stable, tps = full
    assert len(P) == len(T) == len(stable)
    assert tps.shape[1] == 2
    assert _line('Ih', 'VI') is None   # no such pair: returns None, does not raise
