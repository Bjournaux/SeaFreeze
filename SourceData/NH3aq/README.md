# NH3aq: source of the SeaFreeze aqueous-ammonia spline

Everything needed to rebuild `Python/seafreeze/splines/NH3_aq_RK_2026/NH3_aq_RK_2026.mat`.

```bash
python SourceData/NH3aq/build_NH3aq_spline.py
```

Run from the repository root. Needs the `seafreeze` package of this repository (no CoolProp).
The script prints the accuracy of the spline against the analytic model and refuses to
write a spline that misses the water chemical potential by more than 1 J/mol.

| File | Content |
|---|---|
| `nh3aq_model.py` | The model: water reference, liquid NH3, Redlich-Kister excess term, `G_solution(P, T, m)`, and the closed-form water chemical potential `muw_exact` |
| `build_NH3aq_spline.py` | Samples `G_solution`, builds the spline with `seafreeze.splinebuilder`, checks it, writes the `.mat` |
| `rk_parameters.json` | Fitted Redlich-Kister parameters, uncertainties and fit statistics |
| `liquidus_data.csv` | The 88 ice liquidus points used for the fit |
| `nh3_liquid_gibbs_CoolProp.csv` | Molar Gibbs energy of pure liquid NH3, 1-900 MPa, 233.4-273 K |

## Model

Per kg of water (n_w mol of H2O, m mol of NH3; x1, x2 the mole fractions of H2O and NH3):

```
G_kgw = G_w(P,T) + m mu_NH3(P,T) + R T [m ln x2 + n_w ln x1] + (n_w + m) g_ex
g_ex  = R T x1 x2 sum_k L_k(P,T) (x1 - x2)^k
L_k   = a_k + b_k (T_ref/T - 1) + v_k (P - P_ref)/1000        T_ref = 240 K, P_ref = 0.1 MPa
```

The spline stores `G = G_kgw / (1 + M_NH3 m)` in J per kg of solution, as lbftd expects.
Fitted parameters: a0 = -5.53, a1 = +1.26, b0 = -3.08, v0 = +1.32 per GPa; all others zero.

- **Water**: SeaFreeze `water1` for T >= 240 K; below, a constant-heat-capacity continuation
  from 240 K (G, S and Cp continuous).
- **Liquid NH3**: CoolProp `HEOS::Ammonia` Gibbs energies exported by L. Delaroque, fitted with a
  smooth least-squares spline (residual 0.0003 J/mol) and continued outside the table with
  constant Cp in temperature and a second-order expansion in pressure. This term does not
  enter the chemical potential of water.

## Provenance

- Liquidus data digitised by L. Delaroque (IPGP) from Choukroun & Grasset (2010, J. Chem. Phys.
  133, 144502) and Leliwa-Kopystynski, Maruyama & Nakajima (2002, Icarus 159, 518-528).
- The fit is a rework of her Redlich-Kister model: one liquid model for all ices, water from
  `water1`, ices from SeaFreeze, each high-pressure point attributed to the stable ice
  (ice II replaces ice III at low temperature, ice VI replaces ice V at high pressure).
  Fit RMS 1.22 K on 88 points (original model on the same points: 4.9 K).

## Open points on the data

- Leliwa compositions (3.34, 4.72, 17.6, 22.9) are read as wt%; to confirm against the paper.
- `Leliwa_x11.1` in the original folder is a byte-identical copy of `Leliwa_x4.72`; not used.
- `amm014_III` (x = 0.026 as given) lies 9-12 K below the pure-water melting curve; not used.
- `amm013_II` has no composition; not used.
- Whether the low-temperature points labelled "ice III" were ice II in the experiments is not
  known; the 20 wt% series is only consistent with ice II, the 10 wt% series cannot tell.

## Limitations (draft text for the paper)

> The liquid model is an empirical Redlich-Kister expansion calibrated exclusively on ice
> liquidus measurements at x(NH3) <= 0.24, pressures up to 970 MPa and temperatures between
> 210 and 270 K. Each measurement above 200 MPa was attributed to the ice polymorph that is
> thermodynamically stable at its pressure and temperature, so that ice II, rather than ice
> III, bounds the liquidus below about 250 K; the phase actually crystallised in the
> experiments is not always known, and metastable ice III would shift those points by a few
> kelvin. The data constrain the activity of water along the liquidus, which the model
> reproduces with a root-mean-square deviation of 1.2 K, but not the excess Gibbs energy over
> the full composition range: the composition parameters are strongly correlated, and the
> activity of ammonia and all properties at higher ammonia content are extrapolations that
> should not be used. The temperature and pressure dependences are each a single term on the
> leading parameter, equivalent to a constant excess enthalpy and a pressure-independent
> excess volume. Both were obtained from the shape of the liquidi, along which pressure and
> temperature are correlated, and not from calorimetric or volumetric data; they are only
> marginally resolved, and the sign of the pressure term is opposite to the volume contraction
> known for ammonia-water mixing. Heats of mixing, densities, compressibilities and sound
> speeds must therefore not be derived from the model. Below 240 K the Gibbs energy of pure
> liquid water is not available from experiment at these pressures; it is defined here by a
> constant-heat-capacity extrapolation of the Bollengier et al. (2019) equation of state from
> a reference temperature of 240 K. This neglects the heat-capacity anomaly of supercooled
> water below about 50 MPa; using the Holten et al. (2014) heat capacity instead shifts the
> calculated liquidus by less than 0.1 K up to 20 wt% NH3. The fitted parameters are tied to
> this convention and must be used with it. Finally, the liquid is treated as a molecular
> mixture: ammonia hydrates, which replace ice on the liquidus at higher concentration and
> whose stability field expands under pressure, and ionic speciation are not included, so the
> model gives the ice liquidus only and may predict ice where a hydrate is the stable solid.
