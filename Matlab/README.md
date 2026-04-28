# SeaFreeze

V1.1.0 (Matlab version)

The SeaFreeze package computes thermodynamic and elastic properties of water, ice polymorphs (Ih, II, III, V, VI, VII and X) and aqueous NaCl solutions in the 0-100 GPa and 220-10000 K range, with the study of icy worlds and their oceans in mind. It is based on Gibbs Local Basis Function (LBF) parametrizations (https://github.com/jmichaelb/LocalBasisFunction) for each phase. The formalism is described in Brown (2018), Journaux et al. (2020), and in the liquid water Gibbs parametrization by Bollengier, Brown, and Shaw (2019).

## What's new in 1.1.0
- **Aqueous NaCl solutions** (`NaClaq`) via the same 3D (P,T,m) spline used by the Python version, including a corrected scatter-input mixing-quantities path that uses a per-row baseline at `m = cutoff` (previously broken with a runtime warning).
- **Selective property computation**: ask for only the properties you need (e.g. just `rho` or `{'G','Cp'}`) to save time.
- **No Curve Fitting Toolbox required** for `SF_getprop`. A single evaluator (`fnGval` over `sp_val`) handles every phase, including the LBF-form spline used for ice VII/X. (`SF_WhichPhase` and `SF_PhaseLines` still call `fnval` and so do still need the toolbox; tracked for a future port.)
- **New outputs**: `Js` (Joule-Thomson coefficient) and `gamma_Gruneisen` (Grüneisen parameter) for every phase; `f`, `m`, `xs`, `xw`, `Vw` (partial molar volume of solvent) for `NaClaq` mixing.
- **Optional `sp.Tc`** (dimensionless temperature, `tau = log(T/Tc)`) and **`sp.mask`** (validity-domain interpolation) supported by `fnGval` for new spline parametrizations.
- **Cross-validation test suite** comparing the MATLAB output against the Python reference implementation: 14 cases including grid + scatter + edge molalities for NaCl, with per-case relative-tolerance overrides for known ice-V drifts and freshness/manifest checks on the reference `.mat`.

## Getting started

### Prerequisites
Tested on MATLAB R2018a and newer. `SF_getprop` is now toolbox-free — it uses an in-tree de-Boor evaluator (`sp_val.m`) for every phase. The phase-stability helper `SF_WhichPhase` and the phase-line plotter `SF_PhaseLines` still call `fnval` from the **Curve Fitting Toolbox** and require it to run; this will be removed in a future release.

### Installing
Add `SeaFreeze/Matlab` to your MATLAB path. The following data files must be on the path (they ship with this folder):
- `SeaFreeze_Gibbs.mat` — all ice and pure-water splines
- `SeaFreeze_Gibbs_VII_NaCl5GPa.mat` — aqueous NaCl 3D spline (shared with the Python version)

## Running SeaFreeze

```matlab
out = SF_getprop(PT, material)           % all supported properties (default)
out = SF_getprop(PT, material, props)    % only the requested properties
```

> **Deprecation note (1.1.0):** the legacy entry point `SeaFreeze(PT, material, ...)` still works and is kept as a thin alias for `SF_getprop`. It emits a one-time warning per MATLAB session. Suppress it with `warning('off','SeaFreeze:deprecated')`. Migrate calls to `SF_getprop` at your convenience; the alias will be removed in a future release.

### Inputs

**`PT`** — pressure–temperature (–molality) coordinates.
- Pure phases: cell `{P,T}` (gridded output) or N×2 array `[P T]` (scatter output).
- `NaClaq`: cell `{P,T,m}` (gridded) or N×3 array `[P T m]` (scatter).
- Units: P in MPa, T in K, m in mol/kg. Points outside the parametrization return `NaN`.

**`material`**
| Name | Description |
|------|-------------|
| `Ih` | Ice Ih (Feistel and Wagner, 2006) |
| `II`, `III`, `V`, `VI` | Ices II–VI (Journaux et al. 2020) |
| `VII_X_French` | Ice VII / ice X (French and Redmer 2015) |
| `water1` | Liquid water — Bollengier et al. 2019 (≤500 K, ≤2300 MPa) |
| `water2` | Liquid water — Brown 2018 (up to 100 GPa) |
| `water_IAPWS95` | IAPWS95 water (Wagner and Pruss, 2002) |
| `NaClaq` | Aqueous NaCl solution (3D LBF) |

**`props`** *(optional)* — a string or cell array of property names. Omit or pass `[]` to compute all supported properties.

### Outputs
`out` is a struct with the requested fields (SI units).

Common to all materials:

| Quantity | Field | Unit |
|----------|:-----:|:----:|
| Gibbs energy | `G` | J/kg |
| Entropy | `S` | J/K/kg |
| Internal energy | `U` | J/kg |
| Enthalpy | `H` | J/kg |
| Helmholtz free energy | `A` | J/kg |
| Density | `rho` | kg/m³ |
| Isobaric heat capacity | `Cp` | J/kg/K |
| Isochoric heat capacity | `Cv` | J/kg/K |
| Isothermal bulk modulus | `Kt` | MPa |
| Pressure derivative of `Kt` | `Kp` | – |
| Isentropic bulk modulus | `Ks` | MPa |
| Thermal expansivity | `alpha` | 1/K |
| Bulk sound speed | `vel` | m/s |
| Joule-Thomson coefficient | `Js` | K/MPa |
| Grüneisen parameter | `gamma_Gruneisen` | – |

Solid phases additionally provide:

| Quantity | Field | Unit |
|----------|:-----:|:----:|
| Shear modulus | `shear` | MPa |
| P-wave velocity | `Vp` | m/s |
| S-wave velocity | `Vs` | m/s |

`NaClaq` additionally provides mixing properties:

| Quantity | Field | Unit |
|----------|:-----:|:----:|
| Molality echo | `m` | mol/kg |
| Solute mole fraction | `xs` | – |
| Solvent mole fraction | `xw` | – |
| kg-of-solution per kg-of-water factor | `f` | – |
| Solute chemical potential | `mus` | J/mol |
| Solvent chemical potential | `muw` | J/mol |
| Partial molar volume of solute | `Vm` | cm³/mol |
| Partial molar volume of solvent | `Vw` | cm³/mol |
| Partial molar heat capacity of solute | `Cpm` | J/mol/K |
| Apparent molar volume | `Va` | cm³/mol |
| Apparent molar heat capacity | `Cpa` | J/mol/K |
| Excess volume | `Vex` | cm³/mol |
| Osmotic coefficient | `phi` | – |
| Water activity | `aw` | – |

## Examples

An executable live script (`Example_SeaFreeze.mlx`) is provided.

### Single point (all properties)
```matlab
PT = [900 255];
out = SF_getprop(PT, 'VI');
```

### Grid input
Ice V every 2 MPa from 400 to 500 MPa and every 0.5 K from 240 to 250 K:
```matlab
PT = {400:2:500, 240:0.5:250};
out = SF_getprop(PT, 'V');
```

### Scatter input
Three points for liquid water:
```matlab
PT = [200 300; 223 300; 225 300];
out = SF_getprop(PT, 'water1');
```

### Selective properties (new in 1.1.0)
Only density for an ice VI grid:
```matlab
out = SF_getprop({400:10:900, 240:5:270}, 'VI', 'rho');
```
G, ρ and Cp only:
```matlab
out = SF_getprop([900 255], 'VI', {'G','rho','Cp'});
```
Requesting shear/Vp/Vs pulls in the minimum dependencies automatically:
```matlab
out = SF_getprop([900 255], 'VI', 'Vp');   % returns only Vp
```

### Aqueous NaCl (new in 1.1.0)
For `NaClaq`, `PT` is extended to `(P, T, m)` with molality `m` in mol/kg.

Single point — osmotic coefficient, water activity and apparent properties at 200 MPa, 280 K, 0.5 mol/kg NaCl:
```matlab
PTm = [200 280 0.5];
out = SF_getprop(PTm, 'NaClaq', {'phi','aw','Vex','Va','Cpa'});
% out.phi   ≈ 0.9042   (osmotic coefficient)
% out.aw    ≈ 0.9838   (water activity)
% out.Vex   ≈ 0.444    (apparent excess volume, cm^3/mol)
% out.Va    ≈ 23.12    (apparent molar volume, cm^3/mol)
% out.Cpa   ≈ 15.37    (apparent molar heat capacity, J/mol/K)
```

Scatter list — three arbitrary (P,T,m) conditions:
```matlab
PTm = [100 298 0.5; 200 323 1.0; 500 373 3.0];
out = SF_getprop(PTm, 'NaClaq', {'rho','Cp','mus','muw','aw'});
% out.rho, out.Cp, out.aw, ...   each [3x1]
```

Grid input — sweep pressure 0.1–500 MPa, temperature 273–400 K, molality 0.1–3 mol/kg, all properties:
```matlab
P = 0.1:100:500;        % MPa
T = 273:25:400;         % K
m = [0.1 0.5 1.0 3.0];  % mol/kg
out = SF_getprop({P,T,m}, 'NaClaq');
% rows -> P, columns -> T, third dim -> m
% e.g. out.rho is size [length(P) length(T) length(m)]
```

Subset for speed (e.g. only density and water activity on the same grid):
```matlab
out = SF_getprop({P,T,m}, 'NaClaq', {'rho','aw'});
```

`NaN` is returned outside the parametrization bounds (≤5000 MPa, 229–501 K, up to ~7 mol/kg).

## Tests

A cross-validation suite compares the MATLAB output against the Python reference implementation on a set of ice, water and NaCl cases.

### 1. Generate the Python reference
From `Matlab/test/`:
```bash
python gen_python_reference.py
```
This writes `reference_SeaFreeze.mat`. Requires `numpy`, `scipy`, `hdf5storage`, and the repository's `Python/` directory resolvable on `sys.path` (the script handles that automatically).

### 2. Run the MATLAB comparison
From MATLAB, with `Matlab/test/` on the path:
```matlab
test_SeaFreeze_vs_python
test_SeaFreeze_vs_python('rtol', 1e-5, 'atol', 1e-8, 'verbose', true)
```
The test reports `[pass]`, `[FAIL]` and `[info]` per property per case, and throws an error if any comparison exceeds tolerance.

### Expected tolerances / known differences
- Pure ices, liquid water (`water1`, `water2`, `water_IAPWS95`), ice VII/X and `NaClaq` (grid + scatter + edge molalities) all match to ~1e-10 (relative) or better.
- **Ice V grid** has 6 documented `[info]` drifts (`S`, `U`, `H`, `Cv`, `Ks`, `alpha`; max relative ~7e-4). These reflect a small difference between MATLAB's `fnval(fnder(…))` and SciPy's `splev` knot-extrapolation behaviour on this particular spline — not a code bug. The test classifies them as informational via per-case `expected_drift` overrides; they appear in the report as `[info]` and do **not** count as failures.
- The test prints a per-case summary (`pass / info / fail`) and a final table. `error()` is raised only on real `[fail]`s.
- The test guards reference freshness: if `reference_SeaFreeze.mat` is older than `gen_python_reference.py`, you get a warning to regenerate; if any expected case is missing, the test errors out with a clear hint.

## Utility functions

### `SeaFreeze_version`
Return the current version string:
```matlab
SeaFreeze_version
% '1.1.0'
```

### `SF_WhichPhase`
Determine which supported phase is stable at a given (P,T). `PT` has the same format as `SeaFreeze`. Output integers: 0 = liquid, 1 = ice Ih, 2 = II, 3 = III, 5 = V, 6 = VI; `NaN` outside all parametrizations.
```matlab
SF_WhichPhase({300,300})   % -> 0 (liquid water)
```

By default the liquid is pure water (`water1`). Pass `'solute','NaCl'` to use aqueous NaCl as the liquid; `PT` then needs the molality axis as `{P,T,m}` or `[P T m]`. The comparison is made on the chemical potential of water (μw_solution vs `G_ice·M_H2O`), so it captures freezing-point depression by salt:
```matlab
% Stability map for a 2 mol/kg NaCl solution from 0-1000 MPa, 240-300 K
out = SF_WhichPhase({0:10:1000, 240:1:300, 2}, 'solute','NaCl');
```

### `SF_PhaseLines`
Compute phase-boundary coordinates (melting or solid-solid) by intersecting Gibbs surfaces. Provides metastable extensions beyond published triple points.
```matlab
SF_PhaseLines('Ih', 'water1', 'plot', 'meta')
```

### `SF_WPD`
Plot the full water phase diagram.
```matlab
SF_WPD
```

## Important remarks

### Water representation
Ice Gibbs parametrizations are optimized for use with `water1` (Bollengier et al. 2019), particularly for phase-equilibrium calculations. Using other water parametrizations will yield incorrect melting curves. `water2` and `water_IAPWS95` are provided for high-pressure extension and comparison only. Use `water1` for the 200–355 K, ≤2300 MPa range.

### Range of validity
Stability prediction is considered valid down to 130 K (ice VI – ice XV transition). The Ih–II transition may be valid down to 73.4 K (Ih – XI). The VII/X representation extends to 1 TPa (1e6 MPa) and 2000 K.

## References
- [Bollengier, Brown and Shaw (2019) J. Chem. Phys. 151, 054501](https://aip.scitation.org/doi/abs/10.1063/1.5097179)
- [Brown (2018) Fluid Phase Equilibria 463, 18-31](https://www.sciencedirect.com/science/article/pii/S0378381218300530)
- [Feistel and Wagner (2006) J. Phys. Chem. Ref. Data 35, 1021-1047](https://aip.scitation.org/doi/abs/10.1063/1.2183324)
- [Journaux et al. (2020) JGR: Planets 125, e2019JE006176](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2019JE006176)
- [Wagner and Pruss (2002) J. Phys. Chem. Ref. Data 31, 387-535](https://aip.scitation.org/doi/abs/10.1063/1.1461829)
- [French and Redmer (2015) Phys. Rev. B 91, 014308](http://link.aps.org/doi/10.1103/PhysRevB.91.014308)

## Authors

* **Baptiste Journaux** — *University of Washington, Earth and Space Sciences, Seattle, USA*
* **J. Michael Brown** — *University of Washington, Earth and Space Sciences, Seattle, USA*
* **Penny Espinoza** — *University of Washington, Earth and Space Sciences, Seattle, USA*
* **Ula Jones** — *University of Washington, Earth and Space Sciences, Seattle, USA*
* **Erica Clinton** — *University of Washington, Earth and Space Sciences, Seattle, USA*
* **Tyler Gordon** — *University of Washington, Department of Astronomy, Seattle, USA*

## License

SeaFreeze is licensed under the GPL-3 License:

Copyright (c) 2019, B. Journaux

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

## Acknowledgments

This work was produced with the financial support of the NASA Postdoctoral Program fellowship, the NASA Solar System Workings Grant 80NSSC17K0775, and the Icy Worlds node of NASA's Astrobiology Institute (08-NAI5-0021).
