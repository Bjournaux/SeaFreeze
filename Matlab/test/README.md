# SeaFreeze MATLAB tests

Test suites that exercise the MATLAB SeaFreeze package — both internal correctness checks and cross-validation against the Python reference.

## Quick run

From the `Matlab/` directory:

```matlab
run('run_all_tests.m')
```

This runs all script-based tests and auto-discovers classdef unittest suites in `test/`.

Individual suites can also be run directly:

```matlab
% Script-based suites
test_general              % 35 tests : smoke + thermodynamic identities
test_input_validation     % 20 tests : SeaFreeze:badInput / unknownMaterial paths
test_selective_props      % 18 tests : props argument behaviour
test_fnGval_vs_1p0        % 193 tests : new fnGval vs frozen v1.0 baseline
test_SF_PhaseLines        % 63 tests : phase-line rewrite vs frozen v1 baseline

% Classdef suites (auto-discovered by runtests)
runtests('test')          % test_SeaFreeze (34 tests) + test_SF_WhichPhase (20 tests)

% Cross-version against Python (requires reference_SeaFreeze.mat — see below)
test_SeaFreeze_vs_python
```

All suites should be green. Total: **383+ pass / 0 fail**.

## Test files

| File | Purpose |
|---|---|
| [test_general.m](test_general.m) | Smoke tests for all 13 materials, grid/scatter shape checks, thermodynamic identities (Cp−Cv, Ks/Kt = Cp/Cv), NaN propagation, WhichPhase, freezing-point depression. |
| [test_input_validation.m](test_input_validation.m) | Bad inputs raise the expected `SeaFreeze:badInput` / `unknownMaterial` / `unknownProperty` errors. Includes deprecated `SeaFreeze()` alias warning check. |
| [test_selective_props.m](test_selective_props.m) | The `props` argument selects exactly the requested fields; shear/Vp/Vs internal dependencies are stripped. Verifies field counts for all-props mode (water: 17, ice: 20, NaClaq: 31). |
| [test_fnGval_vs_1p0.m](test_fnGval_vs_1p0.m) | Compares the current `fnGval.m` against the frozen `legacy/fnGval_1p0.m` baseline across ices/water1/NaClaq, grid + scatter, at `rtol=1e-10`. Documented `muw`/`aw` overrides at `1e-6` for the deliberate water-molar-mass refinement. |
| [test_SF_PhaseLines.m](test_SF_PhaseLines.m) | Smoke for all supported pairs (pure + NaClaq), v1 regression for the 11 pure pairs, symmetry, NaClaq physics (FPD ≈ 3.31 K at m=1, P=0), multi-molality, error paths, plot path. |
| [test_SeaFreeze.m](test_SeaFreeze.m) | Classdef unittest — VI/III/VII shear/Vp/Vs checks against hardcoded expected values. |
| [test_SF_WhichPhase.m](test_SF_WhichPhase.m) | Classdef unittest — phase identification: single/multi-point scatter, grid, NaCl modes. |
| [test_SeaFreeze_vs_python.m](test_SeaFreeze_vs_python.m) | Cross-version comparison against a Python-generated `reference_SeaFreeze.mat` (see below). Per-case rtol overrides for known ice-V drifts. |

## Utility files

| File | Purpose |
|---|---|
| [gen_getProp_reference.m](gen_getProp_reference.m) | Generates `reference_getProp.mat` for Python cross-validation (`test_getProp_vs_matlab.py`). |
| [gen_phaselines_reference.m](gen_phaselines_reference.m) | Generates `reference_phaselines.mat` for Python cross-validation (`test_phaselines_vs_matlab.py`). |
| [gen_python_reference.py](gen_python_reference.py) | Python script that produces `reference_SeaFreeze.mat` for `test_SeaFreeze_vs_python`. |
| [compare_SF_PhaseLines.m](compare_SF_PhaseLines.m) | Visual diagnostic — produces side-by-side v1 vs new comparison figures (3 PNGs). Not a pass/fail test. |
| [strip_NaClaq_splines.m](strip_NaClaq_splines.m) | One-time script that strips fitting metadata from NaClaq spline .mat files. |
| [timing_benchmark.m](timing_benchmark.m) | Measures spline load times and SF_getprop call times. |

## Cross-validation against Python

Two cross-validation approaches exist:

### 1. Matlab → Python (recommended)

MATLAB generates reference `.mat` files that the Python test suite loads and compares against:

```matlab
% From Matlab/ directory
gen_getProp_reference      % writes test/reference_getProp.mat
gen_phaselines_reference   % writes test/reference_phaselines.mat
```

Then from `Python/`:
```bash
python -m pytest seafreeze/test/test_getProp_vs_matlab.py -v    # 9 cases, 218 subtests
python -m pytest seafreeze/test/test_phaselines_vs_matlab.py -v  # 18 cases
```

### 2. Python → Matlab

A Python-generated reference is compared by MATLAB:

```bash
cd Matlab/test && python gen_python_reference.py   # writes reference_SeaFreeze.mat
```

```matlab
test_SeaFreeze_vs_python
test_SeaFreeze_vs_python('rtol', 1e-5, 'atol', 1e-8, 'verbose', true)
```

### Expected tolerances / known differences

- Pure ices, liquid water, ice VII/X, and NaClaq all match Python to ~1e-10 (relative) or better.
- **Ice V grid** has documented `[info]` drifts (`S`, `U`, `H`, `Cv`, `Ks`, `alpha`; max relative ~7e-4). These reflect a small difference between MATLAB's `fnval(fnder(…))` and SciPy's `splev` knot-extrapolation behaviour — not a code bug.
- NaClaq stitched cases: `Js` and `gamma_Gruneisen` in the blend zone differ by up to ~1.3e-3 due to different spline evaluation backends (lbftd vs fnGval). `Vw` differs by up to ~1.6e-3 (finite-difference vs analytical derivative).

## Phase-line comparison figures

`compare_SF_PhaseLines` is a visual diagnostic, not a pass/fail test:

```matlab
compare_SF_PhaseLines              % saves to current dir
compare_SF_PhaseLines('/tmp')      % saves to /tmp
```

## Directory structure

```
test/
├── test_general.m                 # script-based tests
├── test_input_validation.m
├── test_selective_props.m
├── test_fnGval_vs_1p0.m
├── test_SF_PhaseLines.m
├── test_SeaFreeze.m               # classdef unittest
├── test_SF_WhichPhase.m           # classdef unittest
├── test_SeaFreeze_vs_python.m     # cross-validation
├── gen_getProp_reference.m        # reference generators
├── gen_phaselines_reference.m
├── gen_python_reference.py
├── reference_getProp.mat          # generated reference data
├── reference_phaselines.mat
├── reference_SeaFreeze.mat
├── compare_SF_PhaseLines.m        # visual diagnostics
├── strip_NaClaq_splines.m        # one-time utilities
└── timing_benchmark.m
```

## How the suites relate

```
test_general / test_input_validation / test_selective_props
        │
        └─ exercise SF_getprop directly

test_fnGval_vs_1p0  ──── compares against frozen legacy/fnGval_1p0.m
test_SF_PhaseLines  ──── compares against frozen legacy/SF_PhaseLines_v1.m

test_SeaFreeze (classdef)  /  test_SF_WhichPhase (classdef)
        │
        └─ in test/; auto-discovered by runtests('test')

test_SeaFreeze_vs_python  ──── compares against Python via reference_SeaFreeze.mat
        │
        └─ requires gen_python_reference.py to have been run first

gen_getProp_reference / gen_phaselines_reference
        │
        └─ generate .mat files consumed by Python test suite
```
