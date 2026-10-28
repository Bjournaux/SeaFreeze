"""Interface-contract tests for getProp, parametrised over every material.

These do not check numbers against references; they check the promises the
API makes regardless of the values:

* selective mode returns exactly the requested properties, with internal
  prerequisites (alpha, Cp, Cv, Kt, muw, rho, Ks) stripped when not asked for;
* selective mode and all-properties mode agree;
* scatter output has shape (n,), grid output (nP, nT) or (nP, nT, nm);
* a scatter point equals the matching grid cell;
* P/T/m echoes reproduce the input.
"""

import numpy as np
import pytest

import seafreeze as sf
from seafreeze.seafreeze import defpath as DEFPATH
from seafreeze.test.golden_cases import GOLDEN_AXES
from seafreeze.test.gen_golden_reference import build_scatter, build_grid

MATERIALS = sorted(sf.phases)
SOLIDS = [m for m in MATERIALS if sf.phases[m].shear_mod_parms]
NACL = [m for m in MATERIALS if m.startswith('NaClaq')]
COORD = {'P', 'T', 'm'}


def _all_props(material):
    out = sf.getProp(build_scatter(GOLDEN_AXES[material]), material)
    return vars(out)


@pytest.mark.parametrize('material', MATERIALS)
def test_all_mode_property_set(material):
    """The all-mode property set is the documented one for the material class."""
    keys = set(_all_props(material))
    base = {'G', 'S', 'U', 'H', 'A', 'rho', 'Cp', 'Cv', 'Kt', 'Kp', 'Ks',
            'alpha', 'vel', 'Js', 'gamma_Gruneisen', 'P', 'T'}
    expected = set(base)
    if material in SOLIDS:
        expected |= {'shear', 'Vp', 'Vs'}
    if material in NACL:
        expected |= {'mus', 'muw', 'Vm', 'Cpm', 'phi', 'aw', 'Va', 'Vex', 'Cpa',
                     'm', 'xs', 'xw', 'f', 'Vw'}
    assert keys == expected, f'{material}: unexpected={keys - expected} missing={expected - keys}'
    # lbftd-internal quantities must never leak
    assert not keys & {'V', 'gam', 'Gex', 'PTM'}


@pytest.mark.parametrize('material', MATERIALS)
def test_selective_returns_exactly_requested(material):
    all_props = _all_props(material)
    PTm = build_scatter(GOLDEN_AXES[material])
    for name in sorted(all_props):
        out = sf.getProp(PTm, material, DEFPATH, name)
        assert set(vars(out)) == {name}, \
            f'{material}: requesting {name!r} returned {sorted(vars(out))}'


@pytest.mark.parametrize('material', MATERIALS)
def test_selective_strips_prerequisites(material):
    """Js needs alpha/rho/Cp; shear needs rho/Ks; Vw needs muw. None may leak."""
    PTm = build_scatter(GOLDEN_AXES[material])
    requests = [('Js',), ('gamma_Gruneisen',), ('Js', 'gamma_Gruneisen')]
    if material in SOLIDS:
        requests += [('shear',), ('Vp', 'Vs')]
    if material in NACL:
        requests += [('Vw',), ('Vw', 'xs')]
    for req in requests:
        out = sf.getProp(PTm, material, DEFPATH, *req)
        assert set(vars(out)) == set(req), f'{material} {req}: got {sorted(vars(out))}'


@pytest.mark.parametrize('material', MATERIALS)
def test_selective_matches_all_mode(material):
    for mode, builder in (('scatter', build_scatter), ('grid', build_grid)):
        PTm = builder(GOLDEN_AXES[material])
        full = vars(sf.getProp(PTm, material))
        # one property at a time, and the whole set at once
        for req in [(k,) for k in sorted(full)] + [tuple(sorted(full))]:
            sel = vars(sf.getProp(PTm, material, DEFPATH, *req))
            for k in req:
                np.testing.assert_allclose(
                    np.asarray(sel[k], float), np.asarray(full[k], float),
                    rtol=1e-12, equal_nan=True,
                    err_msg=f'{material}/{mode}: selective {k} != all-mode {k}')


@pytest.mark.parametrize('material', MATERIALS)
def test_output_shapes(material):
    axes = GOLDEN_AXES[material]
    n = len(axes['P'])
    sc = vars(sf.getProp(build_scatter(axes), material))
    for k, v in sc.items():
        assert np.asarray(v).shape == (n,), f'{material}/scatter/{k}: {np.asarray(v).shape}'
    gr = vars(sf.getProp(build_grid(axes), material))
    nT = len(axes['T'])
    full = (n, nT) if axes['m'] is None else (n, nT, len(axes['m']))
    for k, v in gr.items():
        shp = np.asarray(v).shape
        if k in COORD or k in ('xs', 'xw', 'f'):
            # coordinate echoes and m-only quantities may be 1-D or broadcastable
            assert np.broadcast_shapes(shp, full) == full, f'{material}/grid/{k}: {shp}'
        else:
            assert shp == full, f'{material}/grid/{k}: {shp} != {full}'


@pytest.mark.parametrize('material', MATERIALS)
def test_scatter_equals_grid_cell(material):
    """Every scatter point (P_i, T_i[, m_i]) equals grid cell [i, i(, i)]."""
    axes = GOLDEN_AXES[material]
    sc = vars(sf.getProp(build_scatter(axes), material))
    gr = vars(sf.getProp(build_grid(axes), material))
    n = len(axes['P']); nT = len(axes['T'])
    full = (n, nT) if axes['m'] is None else (n, nT, len(axes['m']))
    for k in sc:
        if k in COORD:
            continue
        g = np.broadcast_to(np.asarray(gr[k], float), full)
        for i in range(n):
            cell = g[i, i, i] if g.ndim == 3 else g[i, i]
            np.testing.assert_allclose(
                np.asarray(sc[k], float)[i], cell, rtol=1e-9, equal_nan=True,
                err_msg=f'{material}/{k}: scatter point {i} != grid cell')


@pytest.mark.parametrize('material', MATERIALS)
def test_coordinate_echoes(material):
    axes = GOLDEN_AXES[material]
    sc = vars(sf.getProp(build_scatter(axes), material))
    np.testing.assert_array_equal(np.asarray(sc['P'], float), axes['P'])
    np.testing.assert_array_equal(np.asarray(sc['T'], float), axes['T'])
    if axes['m'] is not None:
        np.testing.assert_array_equal(np.asarray(sc['m'], float), axes['m'])
    gr = vars(sf.getProp(build_grid(axes), material))
    np.testing.assert_array_equal(np.ravel(np.asarray(gr['P'], float)), axes['P'])
    np.testing.assert_array_equal(np.ravel(np.asarray(gr['T'], float)), axes['T'])
    if axes['m'] is not None:
        np.testing.assert_array_equal(np.ravel(np.asarray(gr['m'], float)), axes['m'])


def test_unknown_material_raises():
    PT = np.empty(1, object); PT[0] = (100.0, 250.0)
    with pytest.raises(ValueError, match='not recognized'):
        sf.getProp(PT, 'ice_IX')


def test_unknown_property_raises():
    PT = np.empty(1, object); PT[0] = (100.0, 250.0)
    with pytest.raises(Exception):
        sf.getProp(PT, 'Ih', DEFPATH, 'not_a_property')


def test_public_api_surface():
    """Names exported from the package must not disappear silently."""
    for name in ('getProp', 'whichphase', 'phasenum2phase', 'phases',
                 'phase_lines', 'phase_range', 'wpd', 'rho2P'):
        assert hasattr(sf, name), f'seafreeze.{name} missing from public API'
    assert set(sf.phases) == set(GOLDEN_AXES), 'materials changed; update golden_cases.py'
