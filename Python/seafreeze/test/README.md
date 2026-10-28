# seafreeze test suite

Run from `Python/`:

```bash
python -m pytest -q
```

GUI tests live in `SeaFreezeGUI/tests` and run from the repository root:

```bash
python -m pytest -q SeaFreezeGUI/tests
```

Both run in CI on every push and pull request (`.github/workflows/tests.yml`),
on Python 3.11–3.13.

## What each file protects

| File | Kind | Purpose |
|---|---|---|
| `test_golden_regression.py` | regression | Every material, every property, scatter + grid, compared to `golden_reference.json` at rtol 1e-8. Any numerical change fails here. |
| `test_getProp_contract.py` | contract | Selective mode returns exactly what was asked; selective == all-mode; shapes; scatter == grid cell; coordinate echoes; public API names. |
| `test_domain_edges.py` | edge cases | Out-of-domain → NaN (never crash, never empty); mixed in/out; NaClaq LP+HP stitch is an exact cosine blend and continuous at both taper edges. |
| `test_physics_identities.py` | invariants | H = G + TS, A = G − P/ρ, Ks/Kt = Cp/Cv, vel² = Ks/ρ, Vp/Vs vs moduli, xs + xw = 1, aw vs φ, μw below pure water. Catch unit/sign errors independent of fitted values. |
| `test_whichphase.py` | behaviour | Phase picks, all liquid EOS agree, NaClaq variants, solute validation, agreement with `phase_lines` at the Ih melting curve. |
| `test_phaselines.py`, `test_rho2P.py` | behaviour | Public API behaviour and error handling (incl. fault injection for rho2P). |
| `test_seafreeze.py` | unit | Low-level helpers and historical spot values. |
| `test_*_vs_matlab.py` | parity | Cross-validation against `Matlab/test/reference_*.mat`. Skip if the file is missing; set `SEAFREEZE_REQUIRE_MATLAB_REF=1` (as CI does) to fail instead. |

## Regenerating the golden reference

Only when a numerical change is intended (new spline file, corrected formula):

```bash
python -m seafreeze.test.gen_golden_reference
```

Commit `golden_reference.json` together with the change that motivated it so
the diff documents what moved. Evaluation points are defined in
`golden_cases.py`; adding a material to `seafreeze.phases` requires adding a
case there (a test enforces this).

## Checking the suite itself (mutation testing)

Before a release, break the code on purpose and confirm the suite goes red,
e.g. change `_STITCH_P_LO`, drop a `1e6` unit factor, swap arguments in the
shear call, or stop stripping unrequested properties. Twelve such mutations
were tried when this suite was written; eleven were caught and the twelfth
(`<=` → `<` on the LP mask at exactly 1001 MPa) is numerically equivalent
because the LP weight is zero there. For a systematic run, `mutmut` on
`seafreeze/seafreeze.py` works with this layout.
