# SeaFreeze: issues to fix on the main branch

Found on 2026-11-19/20 while adding the experimental `NH3aq_RK` material (branch `nh3aq-binary`).
None of these is specific to NH3. This file lives in the untracked `NH3 Delaroque/` folder, so it
stays visible whichever branch is checked out.

Status key: **open** = present on `master` and on the branch; **fixed on branch** = corrected on
`nh3aq-binary` (commit `01e5500` or later) but still present on `master`.

## 1. Spline evaluator returns wrong first derivatives for order-8 splines — open (most important)

- **Where:** `Python/mlbspline/eval.py`, `evalMultivarSpline`, which calls `scipy.interpolate.splev`.
- **What:** for a spline of order 8 (degree 7) the *first* derivative is wrong. Values and second
  derivatives are correct. Checked against `scipy.interpolate.BSpline` and against finite
  differences of G; both agree with each other and disagree with `splev`. Probable cause: `splev`
  wraps FITPACK, which supports degree 5 at most.
- **Who is affected:** ice V and ice VI, whose temperature dimension is order 8 (`order [6 8]`).
  Every other shipped spline is order 6 or lower and is not affected.
- **Which properties:** entropy `S` (= −dG/dT) and thermal expansivity `alpha` (uses d²G/dPdT), and
  everything built on them: `H`, `U`, `A`, `Cv`, `Ks`, `vel`, `Js`, `gamma_Gruneisen`, `Vp`.
  Not affected: `G`, `rho`, `Cp`, `Kt`, hence melting curves and phase boundaries.
- **Size (entropy, relative error):**
  - ice V at 600 MPa: up to 6e-5 between 190 and 280 K; 4 % at 290 K (last knot interval 276–300 K).
  - ice VI at 1500 MPa: up to 6e-4 between 200 and 320 K; 0.5 % at 370 K, 5 % at 380 K, 25 % at
    390 K (last knot interval 347–400 K).
- **Fix:** evaluate with `BSpline` (one pass per dimension, as in `NH3 Delaroque/RK_revised/spline_eval.py`).
- **After fixing:** regenerate `golden_reference.json` (only ice V and VI should change), rerun the
  Matlab-parity tests (`test_getProp_vs_matlab.py`), and check whether the Matlab evaluator
  (`Matlab/internal/fnGval.m`, `sp_val.m`) has the same problem — not checked.
- **Side effect seen:** the NH3-H2O test spline `Python/mlbspline/test/NH3_H2O_LBF.mat` is order 8
  in all three dimensions; with the current evaluator its density shows a false jump at each
  pressure knot (e.g. 1068 → 1043 kg/m3 at 89.6 MPa).

## 2. Water activity assumed two ions per formula unit — fixed on branch

- **Where:** `Python/lbftd/statevars.py`, `evalWaterActivity`: `ln(aw) = -2 m phi MWv`.
- **What:** correct only for NaCl. For a solute with one particle per formula unit it returns the
  square of the right value. The Matlab code (`fnGval.m`) already uses `nu`.
- **Fix on branch:** uses the solute's `nu`, cast to float (as `uint8`, `-nu` wraps to 254).
  No change for NaCl.

## 3. Chemistry functions dispatch on the name `NaClaq` — fixed on branch

- **Where:** `seafreeze.py`, `phaselines.py`, `rho2P.py` (`phase.startswith('NaClaq')` and
  `== 'NaClaq'`).
- **Fix on branch:** a material registry (`is_solution`, `register_solution`, …).
- **Still open:** the same hard-coding in Matlab (`sf_material_defs.m`, `sf_load_spline.m`,
  `SF_getprop.m`, `SF_WhichPhase.m`, `SF_PhaseLines.m`, `SF_WPD.m`, `SF_rho2P.m`,
  `SF_phase_range.m`) and in the GUI (`SeaFreezeGUI/core/constants.py`, `core/compute.py`, `app.py`).

## 4. `whichphase(…, solute='NaCl')` documented but not accepted — fixed on branch

- `Python/README.md` says `solute='NaCl'` works; on `master` only material names (`'NaClaq'`, …)
  are accepted. On the branch a solute name is resolved to its default material.

## 5. `whichphase` on a grid accepts only one molality — open

- A grid `[P, T, m]` with more than one molality fails with a shape error
  ("could not broadcast input array from shape (2,2,2) into shape (2,2)"). Scatter input works.

## 6. Water activity and osmotic coefficient are referenced to the dilute cutoff, not to pure water — open

- **Where:** `lbftd/statevars.py`, `evalOsmoticCoeff` uses G at the cutoff molality (0.0002 mol/kg
  for NaCl) as if it were pure water.
- **Size:** ln(aw) is biased by about 5e-5 for a one-particle solute at that cutoff (measured on a
  synthetic ideal solution); expect roughly twice that for NaCl. `muw`, and therefore melting
  curves, are not affected.
- **Options:** evaluate the reference at m = 0, or use a much smaller cutoff (NH3aq_RK uses 1e-6).

## 7. Two different molar masses for NaCl — open

- Registry (`seafreeze.py`): `mNaCl_kgmol = 58.44e-3`. Spline files: `MW = 0.058443`.
- Effect: about 1e-5 on `f` (kg solution per kg water) and on anything using the registry value.

## 8. Low-precision gas constant in lbftd — open

- `lbftd/statevars.py`: `R = 8.3144`; elsewhere 8.314462618. Relative effect 7.5e-6 on `phi`, `aw`.

## 9. `Vex` of NaClaq behaves badly in the dilute limit and across the stitch — open, not investigated

- `getProp(…, 'NaClaq').Vex` diverges at m = 0 and reached about −200 cm3/mol at 20 wt%, 260 K,
  900 MPa when evaluated point by point.
- Computing the apparent molar volume from densities, a 0.01 mol/kg reference shifts by about
  17 cm3/mol above the 500 MPa LP/HP stitch, i.e. the m → 0 limit of the HP spline and the m = 0
  values differ slightly. Could be a real inconsistency at the stitch or only amplification of
  tiny density differences; needs a proper look.

## 10. Out-of-domain grid evaluation returns object-dtype arrays — open, seen once

- `getProp` on a grid entirely outside the spline domain returned arrays that `np.isnan` could
  not handle without an explicit `float` cast (seen with `NH3aq_RK` at 5000 MPa requesting
  `rho`, `muw`). Not checked for NaCl or the pure phases.

## 11. Tests that assumed NaCl — fixed on branch

- `test_physics_identities.py` hard-coded two ions in the `xs` and `aw` identities; they now read
  `nu` from the registry. `test_getProp_contract.py` and `test_physics_identities.py` selected
  solutions with `startswith('NaClaq')`; they now use `is_solution`.

## Not SeaFreeze bugs, but worth remembering

- **`NH3_H2O_LBF.mat` ("V1")** is a full NH3-H2O Gibbs spline (0–800 MPa, 239–401 K, 0–24.5 mol/kg,
  temperature axis `ln(T/450 K)`) used only as a loader test. It matches handbook densities at
  1 bar and gives sound speeds; it becomes unphysical (negative heat capacity) at low temperature,
  low concentration and high pressure. It could supply the density and sound-speed behaviour that
  `NH3aq_RK` lacks (hybrid not started).
- **`NH3aq_RK` open items** are listed in `SourceData/NH3aq/README.md` on the branch (Delaroque's
  agreement before publishing, Matlab and GUI support, unresolved Choukroun & Grasset compositions).
