# SeaFreeze MATLAB tests

Test suites that exercise the MATLAB SeaFreeze package — both internal correctness checks and cross-validation against the Python reference.

## Quick run

From the `Matlab/` directory:

```matlab
run('run_all_tests.m')
```

This runs all script-based tests and auto-discovers classdef unittest suites in `test/`. Under GNU Octave (which has no `matlab.unittest`) the classdef suites are skipped and `sf_verify_octave` runs in their place — see [Running under Octave](#running-under-octave).

Individual suites can also be run directly:

```matlab
% Script-based suites
test_general              % 35 tests : smoke + thermodynamic identities
test_input_validation     % 24 tests : SeaFreeze:badInput / unknownMaterial paths, renamed-material aliases
test_selective_props      % 18 tests : props argument behaviour
test_fnGval_vs_1p0        % 193 tests : new fnGval vs frozen v1.0 baseline
test_SF_PhaseLines        % 63 tests : phase-line rewrite vs frozen v1 baseline
test_SF_rho2P             % 17 tests : pressure-from-density inversion
test_fnFval               % 48 tests : water_Brown2026 Helmholtz fluid (1.2)
test_phase_diagrams       %  5 tests : SF_WPD_PT / SF_WPD_rhoT, triple points (1.2)
test_water3_internals     % 13 tests : psi_val grid mode, spline cache, warnings, SF_WPD with water_Brown2026 (1.2)
test_SeaFreeze_vs_python  % 18 cases, 263 comparisons vs the Python implementation

% Classdef suites (auto-discovered by runtests)
runtests('test')          % test_SeaFreeze (34 tests) + test_SF_WhichPhase (20 tests)

```

All suites should be green: 11 suites, **470 checks + 263 cross-version comparisons, 0 fail** (MATLAB R2025b). Under GNU Octave: 223 checks + 263 comparisons + `sf_verify_octave`; `test_fnGval_vs_1p0` is skipped (its frozen baseline needs the Curve Fitting Toolbox).
`test_SeaFreeze_vs_python` is part of `run_all_tests` since 1.2.

## Test files

| File | Purpose |
|---|---|
| [test_general.m](test_general.m) | Smoke tests for all 13 materials, grid/scatter shape checks, thermodynamic identities (Cp−Cv, Ks/Kt = Cp/Cv), NaN propagation, WhichPhase, freezing-point depression. |
| [test_input_validation.m](test_input_validation.m) | Bad inputs raise the expected `SeaFreeze:badInput` / `unknownMaterial` / `unknownProperty` errors. Includes deprecated `SeaFreeze()` alias warning check. |
| [test_selective_props.m](test_selective_props.m) | The `props` argument selects exactly the requested fields; shear/Vp/Vs internal dependencies are stripped. Verifies field counts for all-props mode (water: 17, ice: 20, NaClaq: 31). |
| [test_fnGval_vs_1p0.m](test_fnGval_vs_1p0.m) | Compares the current `fnGval.m` against the frozen `legacy/fnGval_1p0.m` baseline across ices/water_Bollengier2019/NaClaq, grid + scatter, at `rtol=1e-10`. Documented `muw`/`aw` overrides at `1e-6` for the deliberate water-molar-mass refinement. |
| [test_SF_PhaseLines.m](test_SF_PhaseLines.m) | Smoke for all supported pairs (pure + NaClaq), v1 regression for the 11 pure pairs, symmetry, NaClaq physics (FPD ≈ 3.31 K at m=1, P=0), multi-molality, error paths, plot path. |
| [test_SeaFreeze.m](test_SeaFreeze.m) | Classdef unittest — VI/III/VII shear/Vp/Vs checks against hardcoded expected values. |
| [test_SF_WhichPhase.m](test_SF_WhichPhase.m) | Classdef unittest — phase identification: single/multi-point scatter, grid, NaCl modes. |
| [test_SF_rho2P.m](test_SF_rho2P.m) | `SF_rho2P` round trips for water, ices and NaClaq; shape, NaN propagation and error paths. |
| [test_fnFval.m](test_fnFval.m) | **1.2** — `water_Brown2026`: `psi_val`/`fnFval` against lbf-thermo's `psiH2O_val` (1e-13; scatter, grid, virial continuation, out-of-box NaN), `Kp` vs dKt/dP, P→ρ→P round trips, thermodynamic identities, stable/liquid/vapour branches (only thermodynamically stable roots near T_c), `'input','rhoT'` for water_Brown2026 and the Gibbs splines (water_Bollengier2019 grid, ice VI incl. shear/Vp, NaClaq), `SF_rho2P` / `SF_phase_range`, the plain F(ρ,T) spline path, `SF_WhichPhase` / `SF_PhaseLines` with water_Brown2026, and `SF_coexistence` (saturation vs IAPWS-95, sublimation vs IAPWS R14-08, dilute extension). |
| [test_phase_diagrams.m](test_phase_diagrams.m) | **1.2** — `sf_phase_map` spot checks, `sf_triple_points` vs the literature (7 points), `SF_WPD_PT` / `SF_WPD_rhoT` render, `sf_melt_T_dq2026`. |
| [test_water3_internals.m](test_water3_internals.m) | **1.2** — `psi_val` grid mode == scattered in every regime, the `sf_load_spline` cache, the `SeaFreeze:longRuntime` and `SeaFreeze:diluteExtension` warnings (where, and only where, they should fire), `SF_WPD('liquid','water_Brown2026')`. |
| [test_SeaFreeze_vs_python.m](test_SeaFreeze_vs_python.m) | Cross-version comparison against a Python-generated `reference_SeaFreeze.mat` (see below): 18 cases incl. water_Brown2026 and (ρ,T) input, rtol 1e-6. Per-case rtol overrides for known ice-V drifts. |

## Utility files

| File | Purpose |
|---|---|
| [gen_getProp_reference.m](gen_getProp_reference.m) | Generates `reference_getProp.mat` for Python cross-validation (`test_getProp_vs_matlab.py`). |
| [gen_phaselines_reference.m](gen_phaselines_reference.m) | Generates `reference_phaselines.mat` for Python cross-validation (`test_phaselines_vs_matlab.py`). |
| [gen_python_reference.py](gen_python_reference.py) | Python script that produces `reference_SeaFreeze.mat` for `test_SeaFreeze_vs_python` (through the public `getProp`, 1.2). |
| [gen_water3_reference.m](gen_water3_reference.m) | Writes `fixtures/psi_reference.mat` (lbf-thermo `psiH2O_val` at scattered/grid/edge states; needs lbf-thermo and the Curve Fitting Toolbox) and `fixtures/water3_getprop_reference.mat` (MATLAB `SF_getprop('water_Brown2026')`, `SF_coexistence` and water_Bollengier2019 (ρ,T) results for the Python parity tests). |
| [gen_helmholtz_fixture.m](gen_helmholtz_fixture.m) | Writes `fixtures/water_F_test.mat`, a plain F(ρ,T) B-spline fitted to IAPWS-95, for the generic Helmholtz path. |
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
python -m pytest seafreeze/test/test_getProp_vs_matlab.py -v    # 13 cases (incl. water_Brown2026, (rho,T)), 286 subtests
python -m pytest seafreeze/test/test_phaselines_vs_matlab.py -v  # 18 cases
python -m pytest seafreeze/test/test_helmholtz.py -v             # water_Brown2026 vs fixtures/*.mat (MATLAB and psiH2O_val)
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

- Pure ices, liquid water, ice VII/X, NaClaq and water_Brown2026 (incl. (ρ,T) input) all match Python within the 1e-6 tolerance; water_Brown2026 to ~1e-9.
- **Ice V grid** has documented `[info]` drifts (`S`, `U`, `H`, `Cv`, `Ks`, `alpha`; max relative ~7e-4). These reflect a small difference between MATLAB's `fnval(fnder(…))` and SciPy's `splev` knot-extrapolation behaviour — not a code bug.
- NaClaq stitched cases: `Js` and `gamma_Gruneisen` in the blend zone differ by up to ~1.3e-3 due to different spline evaluation backends (lbftd vs fnGval). `Vw` differs by up to ~1.6e-3 (finite-difference vs analytical derivative).

## Running under Octave

The same test scripts run under GNU Octave:

```bash
cd Matlab && octave --no-gui --eval "run_all_tests"
```

`run_all_tests` detects Octave (`sf_is_octave`) and replaces the classdef suites, which need `matlab.unittest`, by `sf_verify_octave`, which checks every material against MATLAB-generated references and exercises the public API. All script-based suites run unchanged; they are several times slower than in MATLAB (the water_Brown2026 suites take minutes).

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
├── test_SF_rho2P.m
├── test_fnFval.m                  # 1.2: water_Brown2026
├── test_phase_diagrams.m          # 1.2: full phase diagrams
├── test_water3_internals.m        # 1.2: internals, warnings
├── test_SeaFreeze_vs_python.m     # cross-validation
├── gen_getProp_reference.m        # reference generators
├── gen_phaselines_reference.m
├── gen_python_reference.py
├── gen_water3_reference.m         # 1.2
├── gen_helmholtz_fixture.m        # 1.2
├── reference_getProp.mat          # generated reference data
├── reference_phaselines.mat
├── reference_SeaFreeze.mat
├── fixtures/                      # 1.2: psi_reference.mat, water3_getprop_reference.mat, water_F_test.mat
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

test_fnFval / test_phase_diagrams / test_water3_internals   (1.2)
        │
        └─ water_Brown2026 via fnFval/psi_val, vs fixtures/psi_reference.mat (lbf-thermo psiH2O_val)

test_SeaFreeze_vs_python  ──── compares against Python via reference_SeaFreeze.mat
        │
        └─ requires gen_python_reference.py to have been run first

gen_getProp_reference / gen_phaselines_reference
        │
        └─ generate .mat files consumed by Python test suite
```
