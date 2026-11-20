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
| `liquidus_data.csv` | The 104 ice liquidus points used for the fit (88 digitised, 16 from the fits of Table III of Leliwa et al.), with the source ice label and the ice used |
| `ice_label_changes.csv` | Every point whose ice in the fit differs from the source label, with the reason and the misfit both ways |
| `nh3_liquid_gibbs_CoolProp.csv` | Molar Gibbs energy of pure liquid NH3, 1-900 MPa, 233.4-273 K |

## Model

Per kg of water (n_w mol of H2O, m mol of NH3; x1, x2 the mole fractions of H2O and NH3):

```
G_kgw = G_w(P,T) + m mu_NH3(P,T) + R T [m ln x2 + n_w ln x1] + (n_w + m) g_ex
g_ex  = R T x1 x2 sum_k L_k(P,T) (x1 - x2)^k
L_k   = a_k + b_k (T_ref/T - 1) + v_k (P - P_ref)/1000        T_ref = 240 K, P_ref = 0.1 MPa
```

The spline stores `G = G_kgw / (1 + M_NH3 m)` in J per kg of solution, as lbftd expects.
Fitted parameters: a0 = -5.34, a1 = +1.14, b0 = -2.69, v0 = +1.41 per GPa; all others zero.

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
  Fit RMS 1.20 K on 104 points in 18 series (original model on the same points: 4.8 K).

## Phase labels changed with respect to the sources (to document in the paper)

The source label of every point is kept in `liquidus_data.csv` (`ice_label`) and never modified;
`ice_fitted` is the ice used. `ice_label_changes.csv` lists the 15 points where they differ, plus
the three 10 wt% points whose source label was kept after being tested (18 rows), with the misfit
both ways. Decisions by B. Journaux, 2026-11-20.

**Choukroun & Grasset 20 wt% "ice III": relabelled ice II (3 points, 217-397 MPa, 215-225 K)**
- Reason: considered an identification error by the experimentalists.
- Evidence: as ice III the series misses by -4.1, -6.3 and -6.9 K, with a trend in pressure; as
  ice II by +1.6, -0.3 and -0.4 K. Ice II has the lower Gibbs energy at these conditions.
  Choukroun & Grasset (2010) themselves conclude that ice III cannot exist above 5-10 wt% NH3, and
  Leliwa-Kopystynski et al. (2002) obtained ice II first when freezing at 250 MPa near 210 K.
- Hard-coded as `RELABEL` in `run_fit.py`.

**Choukroun & Grasset 10 wt% "ice III": best-fitting ice per point (4 points, 239-390 MPa, 235-243 K)**
- Rule: each point takes, between ice II and ice III, the ice giving the smaller
  |T(model) - T(observed)|. Labels and fit depend on each other, so the fit is repeated until the
  labels stop changing (two passes).
- Result: the three points at 239, 290 and 342 MPa keep the source label, ice III (misfit +2.2,
  +1.0, +0.1 K); the point at 390 MPa becomes ice II (+0.7 K, against -0.9 K as ice III).
- Caveats to state in the paper: the labels are selected on the outcome, so the misfit of these four
  points is a lower bound; and ice II has the lower Gibbs energy at all four, so the three ice III
  points are metastable ice III in the model.
- Hard-coded as `BEST_LABEL` in `run_fit.py`.
- Effect on the fit, compared with treating all four as ice II: RMS 1.23 -> 1.20 K, parameters
  moved by a few percent of their uncertainty (a0 -5.26 -> -5.34, a1 +1.10 -> +1.14,
  b0 -2.74 -> -2.69, v0 unchanged at +1.41).

**Changed by the stable-ice rule (11 points)**
- Points labelled ice III or V and not covered above are fitted with the ice of lowest Gibbs energy
  among II, III, V, VI at the measured pressure and temperature: 7 points V -> VI (Choukroun
  20 wt% above 728 MPa, `amm015` above 674 MPa), 1 point III -> V (`amm015`, 379 MPa), and the
  3 unused `amm014` points.
- Ice Ih labels are kept even above 210 MPa (metastable ice I, as stated by Leliwa et al.).

## Status of the data questions

Settled from Leliwa-Kopystynski et al. (2002), full paper
- **Compositions are mass fractions: 3.34, 4.72, 11.1, 17.6 and 22.9 wt%** (text and Table III;
  the abstract rounds the first to 0.034). The fit reads them correctly.
- **The digitised series are faithful.** Against the paper's own polynomial fits (Table III) the
  3.34, 4.72, 17.6 and 22.9 wt% files agree to 0.03-0.17 K on average (0.4 K at worst).
- **`Leliwa_x11.1` is the 4.72 wt% data** (it matches the 4.72 fit to 0.17 K and misses the
  11.1 fit by 12 K). The paper's 11.1 wt% series was never digitised.
- **They observed ice III at 3.34 and 4.72 wt% only** (130-300 MPa), identified visually: ice III
  sinks, ice I floats. Above 220 MPa their ice I points are metastable, as assumed in the fit.
- **Their compositions are not independent measurements.** 20-50 % of the ammonia evaporated
  during loading, and the concentration of each sample was then inferred from its liquidus
  temperature at atmospheric pressure using the Kargel (1992) diagram. The 1-bar end of every
  Leliwa series therefore lies on that diagram by construction.
- Accuracy stated: temperature 0.1 K; liquidus pressure 1 MPa at low concentration, 15 MPa at
  high concentration.

Settled from the Choukroun & Grasset (2010) abstract
- Compositions are in wt% (model valid for 0-33 wt% NH3).
- Their model finds that "ice III cannot exist at concentrations above 5-10 wt% NH3 depending
  on pressure". Here ice III leaves the liquidus between 8 and 11 wt%.

Added to the fit from Table III of Leliwa et al. (T in K, p in MPa), as values of the published fits, not raw points
- Ice I, 11.1 wt%: `T = 256.1 - 0.0759 p - 0.000145 p^2`.
- Ice III, 3.34 wt%: `T = 235.9 + 0.0755 p - 0.000073 p^2`; 4.72 wt%: `T = 235.3 + 0.0596 p - 0.000036 p^2`
  (130-300 MPa; the table prints all coefficients as positive, the linear term must be positive
  for ice III).
- Before they were added, the model predicted them out of sample to +1.0 K (ice I, 11.1 wt%,
  20-240 MPa), -1.5 K (ice III, 3.34 wt%) and -1.2 K (ice III, 4.72 wt%). Adding them left the
  selected model unchanged and moved every parameter by a small fraction of its uncertainty;
  the residuals after refitting are +0.8, -1.5 and -1.3 K. The paper's pure-water ice III fit
  is itself 0.6 K below SeaFreeze.
- Each added series has as many points as a digitised one (8, 4 and 4), so the fits do not
  outweigh the raw data.

Open
- `amm014_III` (x = 0.026 as given) lies 9-12 K below the pure-water melting curve; not used.
  Its composition must come from the Choukroun & Grasset paper.
- `amm013_II` has no composition; not used. Same source.
- `amm015` is used with x = 0.0741 as in the original code; the write-up says 7.3 %.
- The Choukroun 20 wt% "ice III" series is relabelled ice II and the 10 wt% series is decided
  point by point (see above); worth confirming with the original authors or L. Delaroque which
  phase was identified, and how.

## Open items before release

- Confirm with L. Delaroque that her digitised data (`liquidus_data.csv`), the NH3 table
  extracted from her export (`nh3_liquid_gibbs_CoolProp.csv`) and this rework of her model can
  be published in the public SeaFreeze repository, and how she should be credited.
- Matlab support: `sf_material_defs.m`, `sf_load_spline.m`, `SF_getprop.m`, `SF_WhichPhase.m`,
  `SF_PhaseLines.m`, `SF_WPD.m`, `SF_rho2P.m`, `SF_phase_range.m` test for `NaClaq` by name;
  copy the spline to `Matlab/splines/`.
- GUI support: `SeaFreezeGUI/core/constants.py`, `core/compute.py`, `app.py`; then remove the
  experimental-material exemption in `SeaFreezeGUI/tests/test_compute.py`.

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
