"""Golden-value regression suite for getProp.

Compares the current all-properties output of every material against the
committed ``golden_reference.json`` (see ``gen_golden_reference.py``).

Any numerical change — a spline file swap, an edited constant, a refactor
that alters a derivative path — fails here.  If the change is intended,
regenerate the reference and commit the JSON alongside it.
"""

import json
import os

import numpy as np
import pytest

import seafreeze as sf
from seafreeze.test.golden_cases import GOLDEN_AXES, GOLDEN_RTOL
from seafreeze.test.gen_golden_reference import build_scatter, build_grid

HERE = os.path.dirname(os.path.abspath(__file__))
REF_PATH = os.path.join(HERE, 'golden_reference.json')


@pytest.fixture(scope='module')
def golden():
    if not os.path.exists(REF_PATH):
        pytest.fail(f'Golden reference missing: {REF_PATH}. '
                    f'Run: python -m seafreeze.test.gen_golden_reference')
    with open(REF_PATH) as fh:
        return json.load(fh)['cases']


def _from_json(v):
    """JSON null -> NaN, nested lists -> float ndarray."""
    return np.array(v, dtype=float) if isinstance(v, list) \
        else np.array(np.nan if v is None else v, dtype=float)


def _compare(expected, actual, material, mode):
    exp_keys = set(expected)
    act_keys = set(actual)
    assert exp_keys == act_keys, (
        f'{material}/{mode}: property set changed. '
        f'missing={sorted(exp_keys - act_keys)} new={sorted(act_keys - exp_keys)}')
    for name in sorted(exp_keys):
        e = _from_json(expected[name])
        a = np.asarray(actual[name], dtype=float)
        assert e.shape == a.shape, f'{material}/{mode}/{name}: shape {a.shape} != {e.shape}'
        np.testing.assert_allclose(
            a, e, rtol=GOLDEN_RTOL, atol=0, equal_nan=True,
            err_msg=f'{material}/{mode}/{name} drifted from golden reference')


@pytest.mark.parametrize('material', sorted(GOLDEN_AXES))
def test_golden_scatter(golden, material):
    assert material in golden, f'{material} not in golden reference; regenerate.'
    out = sf.getProp(build_scatter(GOLDEN_AXES[material]), material)
    _compare(golden[material]['scatter'], vars(out), material, 'scatter')


@pytest.mark.parametrize('material', sorted(GOLDEN_AXES))
def test_golden_grid(golden, material):
    assert material in golden, f'{material} not in golden reference; regenerate.'
    out = sf.getProp(build_grid(GOLDEN_AXES[material]), material)
    _compare(golden[material]['grid'], vars(out), material, 'grid')


def test_golden_covers_every_material(golden):
    """Adding a material to seafreeze.phases must add a golden case."""
    assert set(sf.phases) == set(GOLDEN_AXES) == set(golden)
