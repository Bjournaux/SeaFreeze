# SeaFreeze

V1.0.2 (Matlab version)

**An update to 1.0 will be available soon, refer to the python version for latest updates and support of aqueous NaCl solutions**

The SeaFreeze package allows to compute the thermodynamic and elastic properties of water and ice polymorphs (Ih, III, V, VI and ice VII/ice X) in the 0-100 GPa and 220 - 10000K range, with the study of icy worlds and their ocean in mind. It is based on the evaluation of Gibbs Local Basis Functions parametrization (https://github.com/jmichaelb/LocalBasisFunction) for each phase. The formalism is described in more details in Brown (2018), Journaux et al. (2019), and in the liquid water Gibbs parametrization by Bollengier, Brown, and Shaw (2019). 

## Getting Started


### Prerequisites

This Matlab version of SeaFreeze has been tested on MATLAB_R2018a. SeaFreeze does not need any other toolboxes to work. Nonetheless, the "Curve Fitting Toolbox" allows to run the code significantly faster when evaluating a list of points rather than a grid.

### Installing

To install SeaFreeze, either change your active directory to SeaFreeze/Matlab, or add that directory to your Matlab path.  In either case, the SeaFreeze_Gibbs.mat file must be in the same directory as the Matlab functions.

## Running SeaFreeze

### Inputs
To run the `SeaFreeze` function you need to provide pressure (MPa) and temperature (K) coordinates and a material input:

```Matlab
out=SeaFreeze(PT,'material')
```

PT is a structure (gridded output) or array (scatter output) containing pressure-temperature points (MPa and Kelvin).

'material' defines which ice or water to use.  Possibilities:

- 'Ih' for ice Ih (Feistel and Wagner, 2006)
- 'II' for ice II (Journaux et al. 2019)
- 'III' for ice III (Journaux et al. 2019)
- 'V' for ice V (Journaux et al. 2019)
- 'VI' for ice VI (Journaux et al. 2019)
- 'VII_X_French' for ice VII and ice X (French and Redmer 2015)
- 'water1' for Bollengier et al. (2019) LBF extending to 500 K and 2300 MPa
- 'water2' for the modified EOS in Brown 2018 extending to 100 GPa
- 'water_IAPWS95' for IAPWS95 water (Wagner and Pruss, 2002)

### Outputs
out is a structure containing all output quantities (SI units):


| Quantity        |  Symbol in SeaFreeze  |  Unit (SI)  |
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
| Shear modulus     |`shear`| MPa |
| P wave velocity     |`Vp`| m/s |
| S wave velocity     |`Vs`| m/s |
| Bulk sound speed     |`vel`| m/s |

 **NaN values returned when out of parametrization boundaries.**




## Example

An executable matlab live script (Example_SeaFreeze.mlx) is provided allowing to run the following examples.

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

## Functions
### `SeaFreeze_version` 
Print the current version of SeaFreeze
```Matlab
SeaFreeze_version

ans =

    '0.9.3'
```

### `SF_WhichPhase`
Function to determine which of the *supported* phases is stable under given pressure and temperature conditions. The `PT` input format is the same as the main `SeaFreeze` function.

The output of the function is an array with an integer indicating the phase number corresponding to the `PT` input.  The phase number 0 means 
liquid water, phase number 1 means ice Ih, phase number 3 means ice III, etc.  

#### Example
Which phase is stable at 300 MPa and 300 K:
```Matlab
SF_WhichPhase({300,300})
ans = 0
```
Liquid water (0) is predicted to be stable.


### `SF_Phaselines`
The `SF_PhaseLines` function calculates phase boundary coordinates (melting or solid-solid) for all phases included in SeaFreeze. It computes the interesection of Gibbs energies surfaces, and is able to provide metastable extentions beyound published tripple points (Journaux et al. 2020).
#### Inputs

 Phaselines requires two phases input in any order, e.g. :
```Matlab
SF_PhaseLines('Ih', 'water1')
```

 #### Options:
```Matlab
'plot','meta'
```

`'plot'` graphically displays the requested phase equilibria coordinates 
`'meta'` provides part of the metastable extensions. 
These options can be used separately or in conjunction with each other.
  #### Outputs

  The output of `SF_Phaselines` is a series of Pressure (MPa) and Temperature (K) coordinates. Additionally, if 'plot' is selected, a figure appears with plotted coordinates.

  #### Example
Ice Ih melting curve with figure and metastable extention
```Matlab
Phaselines('Ih', 'water1','plot', 'meta')
```

### `SF_WPD`
Plot the full water phase diagram
```Matlab
SF_WPD
```

## Important remarks 
### Water representation
The ices Gibbs parametrizations are optimized to be used with 'water1' Gibbs LBF from Bollengier et al. (2019), specially for phase equilibrium calculation. Using other water parametrization wil lead to incorect melting curves. 'water2' (Brown 2018) and 'water_IAPWS95' (IAPWS95) parametrization are provided for HP extention (up to 100 GPa) and comparison only. The authors recommend the use of 'water1' (Bollengier et al. 2019) for any application in the 200-355 K range and up to 2300 MPa.

### Range of validity
SeaFreeze stability prediction is currently considered valid down to 130K, which correspond to the ice VI - ice XV transition. The ice Ih - II transition is potentially valid down to 73.4 K (ice Ih - ice XI transition). The ice VII and ice X representation extend to 1TPa (1e6 MPa) and 2000K.

## References
- [Bollengier, Brown and Shaw (2019) J. Chem. Phys. 151, 054501; doi: 10.1063/1.5097179](https://aip.scitation.org/doi/abs/10.1063/1.5097179)
- [Brown (2018) Fluid Phase Equilibria 463, pp. 18-31](https://www.sciencedirect.com/science/article/pii/S0378381218300530)
- [Feistel and Wagner (2006), J. Phys. Chem. Ref. Data 35, pp. 1021-1047](https://aip.scitation.org/doi/abs/10.1063/1.2183324)
- [Journaux et al., (2019), in review (available on ArXiv)](https://arxiv.org/abs/1907.09598)
- [Wagner and Pruss (2002), J. Phys. Chem. Ref. Data 31, pp. 387-535](https://aip.scitation.org/doi/abs/10.1063/1.1461829)
- [French and Redmer (2015), Physical Review B 91, 014308](http://link.aps.org/doi/10.1103/PhysRevB.91.014308)

## Authors

* **Baptiste Journaux** - *University of Washington, Earth and Space Sciences Department, Seattle, USA* 
* **J. Michael Brown** - *University of Washington, Earth and Space Sciences Department, Seattle, USA* 
* **Penny Espinoza** - *University of Washington, Earth and Space Sciences Department, Seattle, USA* 
* **Erica Clinton** - *University of Washington, Earth and Space Sciences Department, Seattle, USA*
* **Tyler Gordon** - *University of Washington, Department of Astronomy, Seattle, USA*



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

This work was produced with the financial support provided by the NASA Postdoctoral Program fellowship, by the NASA Solar System Workings Grant 80NSSC17K0775 and by the Icy Worlds node of NASA's Astrobiology Institute (08-NAI5-0021).


