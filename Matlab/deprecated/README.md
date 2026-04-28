# Deprecated MATLAB files

Files in this folder are no longer part of the active SeaFreeze MATLAB package. They are kept for **historical reference and for diff/blame purposes**, but should not be added to the MATLAB path. Nothing in `Matlab/`, `Matlab/test/`, or `Matlab/LocalBasisFunction/` references them.

| File | Replaced by | Notes |
|------|-------------|-------|
| `fnGval20022024.m` | `Matlab/fnGval.m` | 878-LOC snapshot dated 2024-02-20. Filename encodes the timestamp; zero callers anywhere in the package. |
| `fnGval2.m` | `Matlab/fnGval.m` | Older experimental variant from `LocalBasisFunction/`. Declares its function as `fnGval` (not `fnGval2` — filename-vs-declaration mismatch), positional-return signature. Zero callers. |

## Files NOT moved here even though they look "deprecated"

A few items in the package look like older versions but are actively used and remain on the main path:

- **`Matlab/fnGval_1p0.m`** — frozen v1.0 baseline of `fnGval`, used by [`test/test_fnGval_vs_1p0.m`](../test/test_fnGval_vs_1p0.m) for regression comparison. Function renamed internally to `fnGval_1p0` so it cannot be called accidentally.
- **`Matlab/SF_PhaseLines_v1.m`** — frozen v1 baseline of `SF_PhaseLines`, used by [`test/test_SF_PhaseLines.m`](../test/test_SF_PhaseLines.m) and [`test/compare_SF_PhaseLines.m`](../test/compare_SF_PhaseLines.m).
- **`Matlab/SeaFreeze.m`** — deprecated *alias* for `SF_getprop`. Still callable so old user scripts keep working; emits a one-time `SeaFreeze:deprecated` warning per session. Exercised by `test/test_input_validation.m`.
- **`Matlab/LocalBasisFunction/fnGval.m`** and **`Matlab/LocalBasisFunction/sp_val.m`** — older versions of the spline evaluators that are shadowed at runtime by the top-level `Matlab/fnGval.m` / `sp_val.m` (path order in tests guarantees the new ones win), but kept in place because [`CalcWaterHug.m`](../LocalBasisFunction/CalcWaterHug.m) uses the old positional-return signature `[rho, vel, G, Cp, ...] = fnGval(...)` rather than the struct-return interface of the current top-level version.
