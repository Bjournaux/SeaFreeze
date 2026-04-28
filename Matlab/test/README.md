# SeaFreeze MATLAB tests

Test suites that exercise the MATLAB SeaFreeze package — both internal correctness checks and a cross-validation against the Python reference.

## Quick run

From MATLAB, with `Matlab/` on the path:

```matlab
% Function-style suites — just call them
test_general              % 30 tests : smoke + thermodynamic identities
test_input_validation     % 18 tests : SeaFreeze:badInput / unknownMaterial paths
test_selective_props      % 18 tests : props argument behaviour
test_fnGval_vs_1p0        % 193 tests : new fnGval vs frozen v1.0 baseline
test_SF_PhaseLines        % 43 tests : phase-line rewrite vs frozen v1 baseline

% classdef-style suites — use runtests
runtests('test_SeaFreeze')         % 34 tests
runtests('test_SF_WhichPhase')     % 17 tests

% Cross-version against Python (requires reference_SeaFreeze.mat — see below)
test_SeaFreeze_vs_python  % 203 pass / 6 info / 0 fail
```

All suites should be green. Total: **356 pass / 6 informational / 0 fail**.

## Test files

| File | Purpose |
|---|---|
| [test_general.m](test_general.m) | Smoke + thermodynamic identities (Cp−Cv = T·α²·Kt/ρ, Ks/Kt = Cp/Cv, NaN propagation, freezing-point depression). |
| [test_input_validation.m](test_input_validation.m) | Bad inputs raise the expected `SeaFreeze:badInput` / `unknownMaterial` / `unknownProperty` errors. |
| [test_selective_props.m](test_selective_props.m) | The `props` argument selects exactly the requested fields; shear/Vp/Vs internal dependencies are stripped. |
| [test_fnGval_vs_1p0.m](test_fnGval_vs_1p0.m) | Compares the current `fnGval.m` against the frozen `fnGval_1p0.m` baseline across ices/water1/NaClaq, grid + scatter, at `rtol=1e-10`. Documented `muw`/`aw` overrides at `1e-6` for the deliberate water-molar-mass refinement (`1000/18.01528` → `1000/18.015268`). |
| [test_SF_PhaseLines.m](test_SF_PhaseLines.m) | Smoke for all 16 supported pairs, v1 regression for the 11 pure pairs, symmetry, NaClaq physics (FPD ≈ 3.31 K at m=1, P=0), error paths, plot path. |
| [test_SeaFreeze_vs_python.m](test_SeaFreeze_vs_python.m) | Cross-version comparison against a Python-generated `reference_SeaFreeze.mat` (see below). Per-case rtol overrides for known ice-V drifts. |
| [compare_SF_PhaseLines.m](compare_SF_PhaseLines.m) | Visual diagnostic — produces side-by-side v1 vs new comparison figures (3 PNGs). Not a pass/fail test; useful when changing the phase-line code. |
| [gen_python_reference.py](gen_python_reference.py) | Python script that produces `reference_SeaFreeze.mat`. Run once before `test_SeaFreeze_vs_python`. |

The classdef-style legacy tests `test_SeaFreeze.m` and `test_SF_WhichPhase.m` live at the package root (`Matlab/test_SeaFreeze.m`, `Matlab/test_SF_WhichPhase.m`) for historical reasons; they are still active and run with `runtests`.

## Cross-validation against Python

`test_SeaFreeze_vs_python` compares MATLAB output against a Python reference on a set of ice, water and NaCl cases.

### 1. Generate the Python reference

From `Matlab/test/`:
```bash
python gen_python_reference.py
```

This writes `reference_SeaFreeze.mat` (14 cases). Requires `numpy`, `scipy`, `hdf5storage`, and the repository's `Python/` directory resolvable on `sys.path` (the script handles that automatically).

The reference includes:
- Pure phases: `Ih_grid`, `Ih_scatter`, `II_grid`, `III_grid`, `V_grid`, `VI_grid`, `VII_X_French_grid`, `water1_grid`, `water1_scatter`, `water_IAPWS95_grid`, `water2_grid`.
- NaClaq: `NaClaq_grid`, `NaClaq_edges_grid` (m at the cutoff and near upper knot), `NaClaq_scatter`.

### 2. Run the MATLAB comparison

From MATLAB, with `Matlab/test/` on the path:
```matlab
test_SeaFreeze_vs_python
test_SeaFreeze_vs_python('rtol', 1e-5, 'atol', 1e-8, 'verbose', true)
```

The test reports `[pass]`, `[FAIL]` and `[info]` per property per case, prints per-case `pass / info / fail` summaries plus a final table, and throws an error only on real `[FAIL]`s.

### Expected tolerances / known differences

- Pure ices, liquid water (`water1`, `water2`, `water_IAPWS95`), ice VII/X and `NaClaq` (grid + scatter + edge molalities) all match Python to ~1e-10 (relative) or better.
- **Ice V grid** has 6 documented `[info]` drifts (`S`, `U`, `H`, `Cv`, `Ks`, `alpha`; max relative ~7e-4). These reflect a small difference between MATLAB's `fnval(fnder(…))` and SciPy's `splev` knot-extrapolation behaviour on this particular spline — not a code bug. The test classifies them as informational via per-case `expected_drift` overrides; they appear in the report as `[info]` and do **not** count as failures.
- The test prints a per-case summary (`pass / info / fail`) and a final table. `error()` is raised only on real `[fail]`s.
- The test guards reference freshness: if `reference_SeaFreeze.mat` is older than `gen_python_reference.py`, you get a warning to regenerate; if any expected case is missing, the test errors out with a clear hint.

## Phase-line comparison figures

`compare_SF_PhaseLines` is a visual diagnostic, not a pass/fail test. It writes three PNG figures for inspection when you change the phase-line code:

```matlab
compare_SF_PhaseLines              % saves to current dir
compare_SF_PhaseLines('/tmp')      % saves to /tmp
```

- **`SF_PhaseLines_v1_vs_new_pure.png`** — 4×3 grid of all 11 pure pairs. v1 (black) vs new full curve (dotted red) vs new stable portion (solid red), with triple-point markers.
- **`SF_PhaseLines_NaCl_melting_curves.png`** — 4 panels (Ih, III, V, VI ↔ NaClaq) showing freezing-point depression at m = 0.1, 0.5, 1, 2, 4 mol/kg vs the pure-water reference.
- **`SF_PhaseLines_full_diagram.png`** — full water phase diagram, v1 underneath in black, new on top in red. Should look like a perfect overlay.

## How the suites relate

```
test_general / test_input_validation / test_selective_props
        │
        ├─ exercise SF_getprop directly
        │
test_fnGval_vs_1p0  ──── compares against frozen Matlab/fnGval_1p0.m
test_SF_PhaseLines  ──── compares against frozen Matlab/SF_PhaseLines_v1.m
                          (also generates the comparison figures via compare_SF_PhaseLines)
        │
test_SeaFreeze (classdef)  /  test_SF_WhichPhase (classdef)
        │
        └─ at Matlab/ root; called via runtests

test_SeaFreeze_vs_python  ──── compares against Python via reference_SeaFreeze.mat
        │
        └─ requires gen_python_reference.py to have been run first
```

## CI hint

A simple driver that runs every suite and exits with a non-zero status on any failure can be built around:

```matlab
function run_all_tests()
    test_general; test_input_validation; test_selective_props;
    test_fnGval_vs_1p0; test_SF_PhaseLines;
    runtests({'test_SeaFreeze','test_SF_WhichPhase'});
    if exist('reference_SeaFreeze.mat','file')
        test_SeaFreeze_vs_python;
    end
end
```
