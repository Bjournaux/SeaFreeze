"""Thermodynamic identity tests.

These hold for any correct Gibbs-energy EOS regardless of the fitted
coefficients, so they catch unit-conversion, sign and factor errors that a
golden-value comparison would simply bake in.  Units: P in MPa, everything
else SI, hence the 1e6 factors.
"""

import numpy as np
import pytest

import seafreeze as sf
from seafreeze.seafreeze import defpath as DEFPATH
from seafreeze.test.golden_cases import GOLDEN_AXES
from seafreeze.test.gen_golden_reference import build_grid

MATERIALS = sorted(sf.phases)
SOLIDS = [m for m in MATERIALS if sf.phases[m].shear_mod_parms]
NACL = [m for m in MATERIALS if m.startswith('NaClaq')]
RTOL = 1e-6


# Stitched NaClaq blends each property separately in the 499.9–1001 MPa taper,
# so nonlinear identities (A = G − P/ρ, Ks/Kt = Cp/Cv, …) only hold approximately
# there.  Exact identities are checked off-taper; the taper deviation is bounded
# by its own test below.
_PHYS_AXES = dict(GOLDEN_AXES)
_PHYS_AXES['NaClaq'] = dict(GOLDEN_AXES['NaClaq'], P=[100.0, 400.0, 3000.0])


def _grid(material):
    axes = _PHYS_AXES[material]
    out = vars(sf.getProp(build_grid(axes), material))
    P = np.asarray(axes['P'], float); T = np.asarray(axes['T'], float)
    if axes['m'] is None:
        P = P[:, None]; T = T[None, :]
    else:
        P = P[:, None, None]; T = T[None, :, None]
    return out, P, T


def _close(a, b, msg):
    a = np.asarray(a, float); b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    assert ok.any(), f'{msg}: no finite points to compare'
    np.testing.assert_allclose(a[ok], b[ok], rtol=RTOL, err_msg=msg)


@pytest.mark.parametrize('material', MATERIALS)
def test_enthalpy_H_equals_G_plus_TS(material):
    o, P, T = _grid(material)
    _close(o['H'], o['G'] + T * o['S'], f'{material}: H != G + T*S')


@pytest.mark.parametrize('material', MATERIALS)
def test_helmholtz_A_equals_G_minus_PV(material):
    o, P, T = _grid(material)
    _close(o['A'], o['G'] - P * 1e6 / o['rho'], f'{material}: A != G - P/rho')


@pytest.mark.parametrize('material', MATERIALS)
def test_internal_energy_U_equals_H_minus_PV(material):
    o, P, T = _grid(material)
    _close(o['U'], o['H'] - P * 1e6 / o['rho'], f'{material}: U != H - P/rho')


@pytest.mark.parametrize('material', MATERIALS)
def test_bulk_modulus_ratio_equals_heat_capacity_ratio(material):
    o, P, T = _grid(material)
    _close(o['Ks'] / o['Kt'], o['Cp'] / o['Cv'], f'{material}: Ks/Kt != Cp/Cv')


@pytest.mark.parametrize('material', MATERIALS)
def test_sound_speed_from_Ks(material):
    o, P, T = _grid(material)
    _close(o['vel'] ** 2, o['Ks'] * 1e6 / o['rho'], f'{material}: vel^2 != Ks/rho')


@pytest.mark.parametrize('material', MATERIALS)
def test_adiabatic_gradient_Js(material):
    o, P, T = _grid(material)
    _close(o['Js'], T * o['alpha'] * 1e6 / (o['rho'] * o['Cp']), f'{material}: Js formula')


@pytest.mark.parametrize('material', MATERIALS)
def test_gruneisen_parameter(material):
    o, P, T = _grid(material)
    _close(o['gamma_Gruneisen'], o['alpha'] * o['Kt'] * 1e6 / (o['rho'] * o['Cv']),
           f'{material}: gamma_Gruneisen formula')


@pytest.mark.parametrize('material', MATERIALS)
def test_heat_capacities_positive_and_Cp_ge_Cv(material):
    o, P, T = _grid(material)
    ok = np.isfinite(o['Cp']) & np.isfinite(o['Cv'])
    assert (np.asarray(o['Cp'])[ok] > 0).all(), f'{material}: Cp <= 0'
    assert (np.asarray(o['Cp'])[ok] >= np.asarray(o['Cv'])[ok] * (1 - 1e-12)).all(), f'{material}: Cp < Cv'


@pytest.mark.parametrize('material', MATERIALS)
def test_density_increases_with_pressure_isothermally(material):
    o, P, T = _grid(material)
    rho = np.asarray(o['rho'], float)
    d = np.diff(rho, axis=0)
    ok = np.isfinite(d)
    assert (d[ok] > 0).all(), f'{material}: rho not monotonic in P'


@pytest.mark.parametrize('material', SOLIDS)
def test_seismic_velocities_consistent_with_moduli(material):
    o, P, T = _grid(material)
    rho, Ks, mu = o['rho'], o['Ks'], o['shear']
    _close(o['Vs'] ** 2, mu * 1e6 / rho, f'{material}: Vs^2 != shear/rho')
    _close(o['Vp'] ** 2, (Ks + 4.0 / 3.0 * mu) * 1e6 / rho, f'{material}: Vp^2 != (Ks+4/3 mu)/rho')
    vp, vs = np.asarray(o['Vp'], float), np.asarray(o['Vs'], float)
    ok = np.isfinite(vp) & np.isfinite(vs)
    assert (vp[ok] > vs[ok]).all(), f'{material}: Vp <= Vs'


@pytest.mark.parametrize('material', NACL)
def test_solution_composition_identities(material):
    o, P, T = _grid(material)
    m = np.asarray(_PHYS_AXES[material]['m'], float)
    xs = np.ravel(np.asarray(o['xs'], float)); xw = np.ravel(np.asarray(o['xw'], float))
    np.testing.assert_allclose(xs + xw, 1.0, rtol=1e-12, err_msg=f'{material}: xs + xw != 1')
    nw = 1.0 / sf.seafreeze.mH2O_kgmol
    np.testing.assert_allclose(xs, 2 * m / (2 * m + nw), rtol=1e-12, err_msg=f'{material}: xs')
    f = np.ravel(np.asarray(o['f'], float))
    np.testing.assert_allclose(f, 1.0 + sf.phases[material].MW * m, rtol=1e-12, err_msg=f'{material}: f')


@pytest.mark.parametrize('material', NACL)
def test_water_activity_and_osmotic_coefficient(material):
    """aw = exp(-phi * nu * m / nw) by definition of the osmotic coefficient,
    and aw < 1 for any real salt concentration."""
    o, P, T = _grid(material)
    m = np.asarray(_PHYS_AXES[material]['m'], float)[None, None, :]
    nw = 1.0 / sf.seafreeze.mH2O_kgmol
    _close(o['aw'], np.exp(-o['phi'] * 2 * m / nw), f'{material}: aw vs phi')
    aw = np.asarray(o['aw'], float)
    assert (aw[np.isfinite(aw)] < 1.0).all(), f'{material}: aw >= 1 with salt present'


@pytest.mark.parametrize('material', NACL)
def test_water_chemical_potential_vs_pure_water(material):
    """Dissolving salt lowers the chemical potential of water: muw(m>0) must be
    below the pure-water Gibbs energy per mole at the same (P, T) wherever
    water1 is valid."""
    axes = _PHYS_AXES[material]
    P = np.asarray(axes['P'], float); T = np.asarray(axes['T'], float)
    o, _, _ = _grid(material)
    PT = np.empty(2, object); PT[0] = P; PT[1] = T
    G_w = np.asarray(sf.getProp(PT, 'water1', DEFPATH, 'G').G, float) * sf.seafreeze.mH2O_kgmol
    muw = np.asarray(o['muw'], float)
    for k in range(muw.shape[2]):
        d = muw[:, :, k] - G_w
        ok = np.isfinite(d)
        if ok.any():
            assert (d[ok] < 0).all(), f'{material}: muw not below pure water G at m index {k}'


@pytest.mark.parametrize('P', [550.0, 750.0, 950.0])
def test_stitch_taper_identity_deviation_is_small(P):
    """Inside the LP/HP taper the per-property blend breaks exact identities;
    the deviation must stay far below the EOS accuracy (0.1 %)."""
    PTm = np.empty(1, object); PTm[0] = (P, 300.0, 2.0)
    o = vars(sf.getProp(PTm, 'NaClaq'))
    rel = abs(o['A'][0] - (o['G'][0] - P * 1e6 / o['rho'][0])) / abs(o['A'][0])
    assert rel < 1e-3, f'A identity deviates by {rel:.2e} at P={P}'
    rel = abs(o['Ks'][0] / o['Kt'][0] - o['Cp'][0] / o['Cv'][0])
    assert rel < 1e-3, f'Ks/Kt vs Cp/Cv deviates by {rel:.2e} at P={P}'
