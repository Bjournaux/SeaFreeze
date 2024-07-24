# SeaFreeze

[![License: GPL-3.0 License](https://img.shields.io/badge/License-GPL3-blue.svg?style=flat-square)](https://opensource.org/license/gpl-3-0/)
[![Twitter Follow](https://img.shields.io/twitter/follow/B_jour.svg?style=flat-square&logo=twitter&label=Follow)](https://twitter.com/B_jour)
[![GitHub Follow](https://img.shields.io/github/followers/Bjournaux.svg?style=flat-square&logo=github&label=Follow)](https://github.com/Bjournaux)

V1.0 
![Logo](https://bjournaux.files.wordpress.com/2019/07/cover.002.png)


The SeaFreeze package allows to compute the thermodynamic and elastic properties of pure water, ice polymorphs (Ih, II, 
III, V VI and ice VII/ice X) up to 100 GPa and 10,000K and aqueous NaCl solution up to 8GPa and 2,000K. It is based on 
the evaluation of Gibbs Local Basis Functions parametrization (https://github.com/jmichaelb/LocalBasisFunction) for each 
phase, constructed to reproduce thermodynamic measurements. The formalism is described in more details in 
[Journaux et al. (2020)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JE006176), and in the liquid water 
Gibbs parametrization by [Bollengier, Brown, and Shaw (2019)](https://aip.scitation.org/doi/abs/10.1063/1.5097179). 
Aqueous NaCl equation of state publication is in preparation.

Currently, the Python version is the most up to date. The Matlab version is still under beta for version 1.0.

Contact: bjournau (at) uw (dot) edu


## Getting Started


### Installing
Report to the README file for each version (Python or Matlab) for installing SeaFreeze.


## Running SeaFreeze

This section provides basic examples on how to run SeaFreeze. It is using pseudocode, so syntax will change 
depending on the version used. 

### Inputs

**For Pure Water and ices**

To run the SeaFreeze function for ices and pure water you need to provide pressure (MPa) and temperature (K) coordinates and a material input:

```
out=SeaFreeze(PT,'material')
```

**For aqueous solutions**

For obtaining properties aqueous solution of a given concentration in mol/kg, you need to provide pressure (MPa), temperature (K) and concentration coordinates (mol/kg):

```
out=SeaFreeze(PTm,'material')
```


**For single properties**

To improve computational efficiency, a list of specified thermodynamic variables can be calculated by specifying them as inputs to the SeaFreeze function:

```
out=SeaFreeze(PT,'material', 'G', 'rho')
```
PT is a structure (gridded output) or array (scatter output) containing pressure-temperature points (MPa and Kelvin).

* 'material' defines which ice, water, or solution to use.  Possibilities:
* 'Ih' for ice Ih (Feistel and Wagner, 2006)
* 'II' for ice II (Journaux et al. 2020)
* 'III' for ice III (Journaux et al. 2020)
* 'V' for ice V (Journaux et al. 2020)
* 'VI' for ice VI (Journaux et al. 2020)
* 'VII_X_French' for ice VII and ice X (French and Redmer 2015)
* 'water1' for Bollengier et al. (2019) LBF extending to 500 K and 2300 MPa
* 'water2' for the modified EOS in Brown 2018 extending to 100 GPa and 10,000 K
* 'water_IAPWS95' for IAPWS95 water (Wagner and Pruss, 2002)
* 'aq_NaCl' for aqueous NaCl from JM Brown and B Journaux et al. (in prep.) 


### Outputs
out is a structure containing all output quantities (SI units):


| Quantity  (PT and PTm)      |  Symbol in SeaFreeze  |  Unit (SI)  |
| --------------- |:---------------------:| :----------:|
| Gibbs Energy           | `G` | J/kg |
| Entropy                | `S` | J/K/kg |
| Internal Energy        | `U` | J/kg |
| Enthalpy               | `H` | J/kg |
| Helmholtz free energy  | `A` | J/kg |
| Density                |`rho`| kg/m^3 |
|Specific heat capacity at constant pressure|`Cp`| J/kg/K |
|Specific heat capacity at constant volume|`Cv`| J/kg/K |
| Isothermal bulk modulus      |`Kt`| MPa |
|Pressure derivative of the Isothermal bulk modulus|`Kp`| - |
| Isoentropic bulk modulus     |`Ks`| MPa |
| Thermal expansivity     |`alpha`| /K |
| Shear modulus (only for solids)    |`shear`| MPa |
| P wave velocity (only for solids)     |`Vp`| m/s |
| S wave velocity (only for solids)     |`Vs`| m/s |
| Bulk sound speed     |`vel`| m/s |


| Quantity  (PTm only)      |  Symbol in SeaFreeze  |  Unit (SI)  |
| --------------- |:---------------------:| :----------:|
| Solute Chemical Potential           | `mus` | J/mol |
| Solvent Chemical Potential                | `muw` | J/mol |
| Partial Molar Volume        | `Vm` | cc/mol |
| Partial Molar Heat Capacity               | `Cpm` | J/kg/K/mol |
| Apparent Heat Capacity  | `Cpa` | J/kg/K/mol |
| Apparent Volume                |`Va`| cc/mol |
|Excess Volume|`Vex`| cc/mol |
|Osmotic Coefficient|`phi`| -|
| Water Activity      |`aw`| - |
|Activity Coefficient|`gam`| - |
| Excess Gibbs Energy     |`Gex`| J/kg |


 **NaN values returned when out of parametrization boundaries.**



## Examples for pure compound (pure water and ices)

An executable Matlab live script (Example_SeaFreeze.mlx) is provided allowing to run the following examples.

### Single point input

Single point for ice VI at 900 MPa and 255 K. This can be used to check returned thermodynamic properties values.
```Matlab
PT = {900,255};
out=SeaFreeze(PT,'VI')
```
Output :
```Matlab
out = 

  struct with fields:

      rho: 1.3561e+03
       Cp: 2.0054e+03
        G: 7.4677e+05
       Cv: 1.8762e+03
      vel: 3.6759e+03
       Kt: 1.7143e+04
       Ks: 1.8323e+04
       Kp: 6.2751
        S: -1.3827e+03
        U: -2.6951e+05
        H: 8.3090e+04
    alpha: 2.0020e-04
       Vp: 4.5490e+03
       Vs: 2.3207e+03
    shear: 7.3033e+03
```

### Grid  input
Grid of points for ice V every 2 MPa from 400 to 500 MPa and every 0.5 K from 220 to 250 K
```Matlab
PT = {400:2:500,240:0.5:250};
out=SeaFreeze(PT,'V')
```
Output :
```Matlab
out = 

  struct with fields:

      rho: [51×21 double]
       Cp: [51×21 double]
        G: [51×21 double]
       Cv: [51×21 double]
      vel: [51×21 double]
       Kt: [51×21 double]
       Ks: [51×21 double]
       Kp: [51×21 double]
        S: [51×21 double]
        U: [51×21 double]
        H: [51×21 double]
    alpha: [51×21 double]
       Vp: [51×21 double]
       Vs: [51×21 double]
    shear: [51×21 double]
```


### List  input
List of 3 points for liquid water at 300K and 200, 223 and 225 MPa 
```Matlab
PT = ([200 300 ; 223 300 ; 225 300 ]);
out=SeaFreeze(PT,'water1')
```

```Matlab
out = 

  struct with fields:

      rho: [3×1 double]
       Cp: [3×1 double]
        G: [3×1 double]
       Cv: [3×1 double]
      vel: [3×1 double]
       Kt: [3×1 double]
       Ks: [3×1 double]
       Kp: [3×1 double]
        S: [3×1 double]
        U: [3×1 double]
        H: [3×1 double]
    alpha: [3×1 double]
```
## Example for Solutions

Thermodynamic properties can be calculated for solutions of varying molality as well, where the input provides pressure (MPa), temperature (K), and molality (mol/kg) coordinates over a grid or list, or for a single point. 

### Single point input

Single point for NaCl(aq) of 0.5 M at 900 MPa and 280 K:

```Matlab
PTm = {900, 280, 0.5};
out=SeaFreeze(PTm,'aq_NaCl')
```
Output :

```Matlab
out = 

  struct with fields:

           G: 7.7725e+05
           S: -143.0425
           U: 1.7738e+04
           H: 7.3720e+05
           A: 5.7790e+04
           F: 5.7790e+04
         rho: 1.2509e+03
          Cp: 3.7630e+03
          Cv: 3.3825e+03
          Kt: 7.7225e+03
          Ks: 8.5913e+03
          Kp: 5.7676
       alpha: 4.6921e-04
         vel: 2.6207e+03
          Va: 24.7939
         Cpa: 124.5329
         mus: 1.7721e+04
         muw: 1.4252e+04
          Vm: 24.3212
          Vw: 14.6032
         Cpm: 149.3719
         gam: 0.7631
         phi: 0.8174
         Vex: 0.2598
         Gex: -204.1910
          aw: 0.9854
```

### Grid input

Grid of points every 10 MPa from 0.1 to 1000 MPa, every 2 K from 240 to 501 K, and every 0.5 M from 1 to 6 mol/kg:

```Matlab
PTm = {0.1:10:1000.2,240:2:501,1:0.5:6}; 
out=SeaFreeze(PTm,'NaCl')
```

Output :

```Matlab

      out = 

  struct with fields:

           G: [101×131×11 double]
           S: [101×131×11 double]
           U: [101×131×11 double]
           H: [101×131×11 double]
           A: [101×131×11 double]
           F: [101×131×11 double]
         rho: [101×131×11 double]
          Cp: [101×131×11 double]
          Cv: [101×131×11 double]
          Kt: [101×131×11 double]
          Ks: [101×131×11 double]
          Kp: [101×131×11 double]
       alpha: [101×131×11 double]
         vel: [101×131×11 double]
          Va: [101×131×11 double]
         Cpa: [101×131×11 double]
         mus: [101×131×11 double]
         muw: [101×131×11 double]
          Vm: [101×131×11 double]
          Vw: [101×131×11 double]
         Cpm: [101×131×11 double]
         gam: [101×131×11 double]
         phi: [101×131×11 double]
         Vex: [101×131×11 double]
         Gex: [101×131×11 double]
          aw: [101×131×11 double]

```

## Important remarks 
### Water representations
The ices' Gibbs parametrizations are optimized to be used with 'water1' Gibbs LBF from Bollengier et al. (2019), 
specially for phase equilibrium calculation. Using other water parametrization wil lead to incorrect melting curves. 
'water2' (Brown 2018) and 'water_IAPWS95' (IAPWS95) parametrization are provided for HP extension (up to 100 GPa) and 
comparison only. The authors recommend the use of 'water1' (Bollengier et al. 2019) for any application in the 200-355 K 
range and up to 2300 MPa.

A Gibbs energy representation of French and Redmer (2015) ice VII and X equation of state is provided for comparison only. It should not be used for melting point or solid-solid phase boundaries predictions.

### Range of validity
SeaFreeze stability prediction is currently considered valid down to 130K, which correspond to the ice VI - ice XV transition. The ice Ih - II transition is potentially valid down to 73.4 K (ice Ih - ice XI transition).




The following figure shows the prediction of phase transitions from SeaFreeze (melting & solid-solid) and comparison with experimental data:
![Logo](https://bjournaux.files.wordpress.com/2019/10/phase-diagram.png)



## Reference to cite to use SeaFreeze:
- [Journaux et al. (2020) JGR Planets 125(1), e2019JE006176 ](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2019JE006176)

## References for the liquid water equations of states used:
- [Bollengier, Brown and Shaw (2019) J. Chem. Phys. 151, 054501; doi: 10.1063/1.5097179](https://aip.scitation.org/doi/abs/10.1063/1.5097179)
- [Brown (2018) Fluid Phase Equilibria 463, pp. 18-31](https://www.sciencedirect.com/science/article/pii/S0378381218300530)
- [Feistel and Wagner (2006), J. Phys. Chem. Ref. Data 35, pp. 1021-1047](https://aip.scitation.org/doi/abs/10.1063/1.2183324)
- [Wagner and Pruss (2002), J. Phys. Chem. Ref. Data 31, pp. 387-535](https://aip.scitation.org/doi/abs/10.1063/1.1461829)
- [French and Redmer (2015), Physical Review B 91, 014308](http://link.aps.org/doi/10.1103/PhysRevB.91.014308)

## Contributors

* **Baptiste Journaux (Lead)** - *University of Washington, Earth and Space Sciences Department, Seattle, USA* 
* **J. Michael Brown** - *University of Washington, Earth and Space Sciences Department, Seattle, USA* 
* **Penny Espinoza** - *University of Washington, Earth and Space Sciences Department, Seattle, USA*
* **Ula Jones** - *University of Washington, Earth and Space Sciences Department, Seattle, USA*
* **Erica Clinton** - *University of Washington, Earth and Space Sciences Department, Seattle, USA*  
* **Tyler Gordon** - *University of Washington, Department of Astronomy, Seattle, USA*
* **Matthew J. Powell-Palm** - *Texas A&M University, Department of Mechanical Engineering, USA*

## Change log

### Changes since 0.9.0
- `1.0`: added NaCl aqueous solution EOS and concentration dependent thermodynamic variables.  NaN returned for values outside the range of the representation.
- `0.9.4`: add ice VII and ice X from French and Redmer (2015).
- [SeaFreeze GUI](https://github.com/Bjournaux/SeaFreeze/tree/master/SeaFreezeGUI) available
- `0.9.4`: Adjusted python readme syntax and package authorship info 
- `0.9.3`: LocalBasisFunction spline interpretation software integrated into SeaFreeze Python package. Adjusted packaging to work better with pip
- `0.9.2` patch1: added `SF_WPD`, `SF_PhaseLines` and `SeaFreeze_version` to the Matlab distribution.
- `0.9.2`: add ice II to the representation.
- `0.9.1`: add `whichphase` function to show which phase is stable at a PT coordinate.

### Planned updates
- Ice VII and X [available here as a beta](https://github.com/Bjournaux/SeaFreeze/tree/new_tdvs) 
- NaCl aqueous solutions [available here as a beta ](https://github.com/Bjournaux/SeaFreeze/tree/new_tdvs/Python) 
- MgSO4, NasSO4 and MgCl2 aqueous solutions
- NH_3 aqueous solutions
- NaCl bearing solids (Halite and hydrohalite)


## License

SeaFreeze is licensed under the GPL-3 License :

Copyright (c) 2019, B. Journaux

This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3.
    
This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY
APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT
HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY
OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM
IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF
ALL NECESSARY SERVICING, REPAIR OR CORRECTION.

## Acknowledgments

As of V1.0, SeaFreeze incorporates the mlbspline and lbftd packages originally developed by J. Michael Brown. Historical versions of these packages are no longer being updated and are available at https://github.com/jmichaelb/LocalBasisFunction. 

This work was produced with the financial support provided by the NASA Postdoctoral Program fellowship, by the NASA Solar System Workings Grant 80NSSC17K0775 and by the Icy Worlds node of NASA's Astrobiology Institute (08-NAI5-0021).

Illustration montage uses pictures from NASA Galileo and Cassini spacecrafts (from top to bottom: Enceladus, Europa and Ganymede). Terrestrial sea ice picture use with the authorization of the author [Rowan Romeyn](https://arcex.no/meet-rowan-romeyn-a-new-arcex-phd-student/).
