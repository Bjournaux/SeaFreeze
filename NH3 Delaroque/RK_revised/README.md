# Revised Redlich-Kister model for the H2O–NH3 ice liquidus

Rework of the model in `../For_Baptiste/` (L. Delaroque). That folder is not modified; its
literature CSVs are read in place. Nothing outside `RK_revised/` is touched.

## Run

```bash
python run_fit.py
```

Needs `numpy`, `scipy`, `pandas`, `matplotlib` and the installed `seafreeze` package.
CoolProp is optional (only for a total liquid Gibbs energy including pure NH3).

| File | Content |
|---|---|
| `nh3_rk.py` | Model: water reference, fusion terms, RK liquid, liquidus solver |
| `run_fit.py` | Data, global fit, model selection, bootstrap, figures, exports |
| `outputs/fit_parameters.json` | Parameters, uncertainties, conventions |
| `outputs/model_ladder.csv` | All candidate models compared |
| `outputs/misfit_by_dataset.csv` | T(model) − T(obs) per series: original, ideal, revised |
| `outputs/data_with_residuals.csv` | Every data point with residuals and ice used |
| `outputs/liquidus_*wt.csv` | Stable liquidus, ice, a_w, ln γ_w on 0–900 MPa |
| `outputs/*.png, *.pdf` | Phase diagram, residuals, water reference |

## Model

One liquid model for all ices (Ih, II, III, V, VI):

```
G_ex = R T x1 x2 Σ_k L_k (x1 − x2)^k                          x1 = x(H2O), x2 = x(NH3)
L_k(P,T) = a_k + b_k (T_ref/T − 1) + v_k (P − P_ref)/1000     T_ref = 240 K, P_ref = 0.1 MPa
ln γ_w = x2² [ L0 + L1 (4x1 − 1) + L2 (12x1² − 8x1 + 1) ]
Liquidus:  G_ice(P,T) − G°_w(P,T) = R T ln(x1 γ_w)
```

**Pure phases**
- Ices: `seafreeze.getProp`, evaluated directly. Outside an EOS domain the result is NaN, never clipped.
- Liquid water, T ≥ 240 K: SeaFreeze `water1` (Bollengier et al. 2019), the liquid the ices are calibrated against.
- Liquid water, T < 240 K: constant-Cp extrapolation of `water1` from 240 K at each pressure,
  `G = G_ref − S_ref (T − T_ref) + Cp_ref [(T − T_ref) − T ln(T/T_ref)]`.
  G, S and Cp are continuous at 240 K (checked). All of this is in `water_ref()` in `nh3_rk.py`,
  so another heat capacity can be swapped in there.
- This choice was tested against the Holten et al. (2014) supercooled-water heat capacity
  (IAPWS G12-15) integrated from the same 240 K anchor: same model selected, RMS 1.224 vs 1.226 K,
  parameters unchanged within a small fraction of their uncertainty, liquidus shifted by
  < 0.1 K up to 20 wt% and by 0.34 K at 22.9 wt% near 1 bar. The Holten option was therefore not kept.

## Result

Selected by leave-one-series-out prediction among 13 parameter subsets: **a0, a1, b0, v0** (all other parameters zero).

| Parameter | Value | Bootstrap 1σ | 95 % interval |
|---|---|---|---|
| a0 | −5.34 | 6.00 | −13.6 to −0.3 |
| a1 | +1.14 | 2.39 | −1.2 to +4.6 |
| b0 | −2.69 | 3.17 | −11.5 to +3.0 |
| v0 (per GPa) | +1.41 | 1.41 | −0.6 to +5.9 |

104 points, 18 series, x = 0.021–0.239, 9–971 MPa, 210.5–270.4 K. Sixteen of the points (three series) are values of the polynomial fits published in Table III of Leliwa-Kopystyński et al. (2002): ice I at 11.1 wt% and ice III at 3.34 and 4.72 wt% (`data_added/`).

| RMS of T(model) − T(obs), fitted points | Original | Revised |
|---|---|---|
| All | 4.76 K | 1.20 K |
| Ice Ih | 4.74 K | 0.91 K |
| Ice II | 2.70 K | 0.91 K |
| Ice III | 3.11 K | 1.39 K |
| Ice V | 4.24 K | 2.06 K |
| Ice VI | 7.99 K | 1.47 K |

"Original" is the `For_Baptiste` code refitted as delivered and evaluated on the same points with the same ice.

Adding the Table III series (they had been predicted out of sample to +1.0, −1.5 and −1.2 K) left the selected model unchanged and moved the parameters by a small fraction of their uncertainty (fit on the original 88 points: a0 −5.53, a1 +1.26, b0 −3.08, v0 +1.32, RMS 1.22 K; on 104 points with all four 10 wt% points as ice II: a0 −5.26, a1 +1.10, b0 −2.74, v0 +1.41, RMS 1.23 K).

**Ice II on the liquidus.** Ice II is the stable ice below the II–III–V triple point (249.4 K, 355 MPa), so it takes over the liquidus from ice III as ammonia lowers the temperature. In the model, ice III survives on the liquidus only up to about 8 wt% (210–255 MPa); from 11 wt% the sequence is Ih → II → V → VI:

| NH3 | Stable liquidus sequence (onset pressure, MPa) |
|---|---|
| 8 wt% | Ih, III (210), II (250), V (380), VI (655) |
| 11 wt% | Ih, II (210), V (405), VI (660) |
| 15 wt% | Ih, II (205), V (445), VI (665) |
| 20 wt% | Ih, II (195), V (510), VI (675) |
| 22.9 wt% | Ih, II (185), V (555), VI (680) |

Fitting the low-temperature "ice III" data as ice II is what the data prefer: the Choukroun 20 wt% series goes from a mean misfit of −3.3 K (as ice III) to +0.3 K (as ice II), and the overall RMS from 1.75 K to 1.22 K. With ice II in place, a temperature term (b0) becomes supported by out-of-sample prediction, and it has the expected sign (exothermic mixing).

Things to know about this result:
- **a0 and a1 are almost perfectly anti-correlated (−1.00)**, and b0 is correlated with both (0.5). Over x ≤ 0.24 the data constrain combinations, not individual values. Use the parameters as a set.
- **b0 and v0 are each within 1–1.1σ of zero.** They improve prediction together, but neither is individually well resolved.
- **The pressure term has an unphysical sign as a volume.** v0 > 0 implies a positive excess volume (+0.24 cm³/mol at x = 0.1), whereas ammonia–water mixing contracts. Treat v0 as an empirical correction to the liquidus slope, not as volumetric information.
- **Ice V is now the worst-fitted phase**; the 7 and 10 wt% series are off in opposite directions (−2.5 K and +1.4 K), which no smooth model can reconcile.
- **Phase labels of the Choukroun "ice III" series are set by decision** (see Data decisions): 20 wt% relabelled ice II, 10 wt% chosen point by point. Recorded in `outputs/ice_label_changes.csv`. Both are independent of `ICE_ASSIGNMENT`; to test the original labels, empty `RELABEL` and `BEST_LABEL` and set `ICE_ASSIGNMENT = "label"`.

## Data decisions

| Decision | Reason |
|---|---|
| `Leliwa_x11.1.csv` not used | Byte-identical copy of `Leliwa_x4.72.csv`. The paper does have an 11.1 wt% series, so it needs re-digitising. |
| Leliwa file numbers read as wt% | **Confirmed**: the abstract of Leliwa-Kopystyński et al. (2002) gives concentrations by mass, X = 0.034, 0.0472, 0.111, 0.176, 0.229. The original code converted only 17.6. (wt% also fits slightly better than mol%: RMS 1.26 vs 1.39 K.) |
| `amm014_III` shown but not fitted | 9–12 K below the pure-water melting curve of ice III or II; impossible for the stated x = 0.026. Composition is wrong. |
| `amm013_II*` not used | No composition given; two inconsistent digitisations |
| Choukroun 20 wt% "ice III" series relabelled ice II (`RELABEL` in `run_fit.py`) | Decision of B. Journaux, 2026-11-20: suspected identification error by the experimentalists. As ice III the series misses by −4 to −7 K, as ice II by +1.6 to −0.4 K. |
| Choukroun 10 wt% "ice III" series: best-fitting ice per point (`BEST_LABEL` in `run_fit.py`) | Decision of B. Journaux, 2026-11-20. Between ice II and III each point takes the ice with the smaller misfit; iterated with the fit until stable. Result: 239, 290, 342 MPa stay ice III; 390 MPa becomes ice II. Labels are selected on the outcome, so the misfit of these points is a lower bound. All logged in `outputs/ice_label_changes.csv`; to document in the paper. |
| Other high-pressure points fitted with the stable ice (`ICE_ASSIGNMENT = "stable"`) | Lowest Gibbs energy among II, III, V, VI at the measured (P, T): 7 points V → VI, 1 point III → V |
| Ice Ih labels kept | The three 2 wt% points at 220–250 MPa follow the metastable Ih melting curve (T still falling with P); fitted as ice III they miss by 7 K |
| `amm015` composition x = 0.0741 kept as in the original | Provenance unknown; the write-up says 7.3 % |
| Added: Leliwa 22.9, Choukroun 10 and 20 wt% Ih, 20 wt% III | Present in the folder but unused in the original |

## Issues found in the original (`For_Baptiste`)

1. Duplicated Leliwa file entered the fit twice with two compositions.
2. L0 ended on its bound (−20), uncertainties far larger than the values; fit still reported success.
3. Mixed composition units for the Leliwa series.
4. Ice V pressure coefficient fitted, then never used; defined per ice although it is a liquid property.
5. IAPWS-95 water tables: negative Cp below ~230 K at ≥ 200 MPa, and 27–45 J/mol away from `water1` at 240 K, 400–600 MPa (see `outputs/water_reference.png`).
6. Table lookups clipped silently outside their P–T range.
7. Exported `G_liquid`: 785 grid points dropped silently and 915 rows with wrong CoolProp water values below 233 K; exported G_ex not constrained by the data.
8. Write-up: Eq. 4 sign (`+8x1` should be `−8x1`; the code was right), Eq. 10 "A, B", eutectic stated at 20 %, Choukroun 2010 reference points to the methane-clathrate paper (should be Choukroun & Grasset 2010, J. Chem. Phys. 133, 144502), SeaFreeze citation.

## Limitations paragraph (draft for the paper)

> The liquid model is an empirical Redlich–Kister expansion calibrated exclusively on ice liquidus measurements at x(NH₃) ≤ 0.24, pressures up to 970 MPa and temperatures between 210 and 270 K. Each measurement above 200 MPa was attributed to the ice polymorph that is thermodynamically stable at its pressure and temperature, so that ice II, rather than ice III, bounds the liquidus below about 250 K; the phase actually crystallised in the experiments is not always known, and metastable ice III would shift those points by a few kelvin. The data constrain the activity of water along the liquidus, which the model reproduces with a root-mean-square deviation of 1.2 K, but not the excess Gibbs energy over the full composition range: the composition parameters are strongly correlated, and the activity of ammonia and all properties at higher ammonia content are extrapolations that should not be used. The temperature and pressure dependences are each a single term on the leading parameter, equivalent to a constant excess enthalpy and a pressure-independent excess volume. Both were obtained from the shape of the liquidi, along which pressure and temperature are correlated, and not from calorimetric or volumetric data; they are only marginally resolved, and the sign of the pressure term is opposite to the volume contraction known for ammonia–water mixing. Heats of mixing, densities, compressibilities and sound speeds must therefore not be derived from the model. Below 240 K the Gibbs energy of pure liquid water is not available from experiment at these pressures; it is defined here by a constant-heat-capacity extrapolation of the Bollengier et al. (2019) equation of state from a reference temperature of 240 K. This neglects the heat-capacity anomaly of supercooled water below about 50 MPa; using the Holten et al. (2014) heat capacity instead shifts the calculated liquidus by less than 0.1 K up to 20 wt% NH₃. The fitted parameters are tied to this convention and must be used with it. Finally, the liquid is treated as a molecular mixture: ammonia hydrates, which replace ice on the liquidus at higher concentration and whose stability field expands under pressure, and ionic speciation are not included, so the model gives the ice liquidus only and may predict ice where a hydrate is the stable solid.

## Not done

- No density, sound-speed or heat-of-mixing data used.
- Holten heat capacity below 240 K tested and removed (no significant effect); NaClaq(m = 0) not tried.
- Total liquid Gibbs energy with pure NH3 not exported (CoolProp not installed here).
- Leliwa compositions are confirmed from the full paper: 3.34, 4.72, 11.1, 17.6, 22.9 wt% (the abstract rounds the first). The digitised files match the paper's Table III fits to 0.4 K or better; `Leliwa_x11.1` is the 4.72 wt% data. Their compositions were inferred from the 1-bar liquidus (Kargel 1992), not measured independently. Table III also gives fits for ice I at 11.1 wt% and ice III at 3.34 and 4.72 wt%, not yet used; the current model predicts them to +1.0 K and −1.2 to −1.5 K. The `amm013`/`amm014`/`amm015` compositions are not in the Choukroun & Grasset abstract and still need the paper's tables.
- Choukroun & Grasset (2010) state that ice III cannot exist above 5–10 wt% NH3, consistent with ice II taking over the liquidus here between 8 and 11 wt%.
