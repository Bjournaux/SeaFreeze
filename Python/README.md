# SeaFreeze

V1.0

The SeaFreeze package allows to compute the thermodynamic and elastic properties of water and ice polymorphs (Ih, III, V, VI and ice VII/ice X) in the 0-100 GPa and 220 - 10000K range, with the study of icy worlds and their ocean in mind. It is based on the evaluation of Gibbs Local Basis Functions parametrization (https://github.com/jmichaelb/LocalBasisFunction) for each phase. The formalism is described in more details in Brown (2018), Journaux et al. (2019), and in the liquid water Gibbs parametrization by Bollengier, Brown, and Shaw (2019). 


## Installation
This package will install SeaFreeze, LBFTD, and MLBspline and their dependencies.

Run the following command to install

`pip install SeaFreeze`

To upgrade to the latest version, use

`pip install --upgrade SeaFreeze`


### `seafreeze.getProps`:
Calculates thermodynamic and elastic properties of a phase of water or solution

### Usage
The main function of SeaFreeze is `seafreeze.getProps`*, which has the following parameters:
- `PT`: the pressure (MPa) and temperature (K) conditions at which the thermodynamic quantities should be
  calculated -- note that these are required units, as conversions are built into several calculations
  This parameter can have one of the following formats:
  - a 1-dimensional numpy array of tuples with one or more scattered (P,T) tuples 
  - a numpy array with 2 nested numpy arrays, the first with pressures and the second
    with temperatures -- each inner array must be sorted from low to high values
    a grid will be constructed from the P and T arrays such that each row of the output
    will correspond to a pressure and each column to a temperature 
- `phase`: indicates the phase of H₂O.  Supported phases are
  - 'Ih' - from  Feistel and Wagner, 2006
  - 'II' - from Journaux et al., 2019
  - 'III' - from Journaux et al., 2019
  - 'V' - from Journaux et al., 2019
  - 'VI' - from Journaux et al., 2019
  - 'VII_X_French' for ice VII and ice X; from French and Redmer 2015
  - 'water1' -  extends to 500 K and 2300 MPa; from Bollengier et al. 2019
  - 'water2' -  extends to 100 GPa; from Brown 2018
  - 'water_IAPWS95' - LBF representation of IAPWS 95; from Wagner and Pruß, 2002
  - 'NaClaq' - LBF representation of aqueous NaCl up to 1000 MPa; from Brown and Journaux in prep

  *Note seafreeze.seafreeze is currently deprecated and will remain so for 2 years as per Python guidelines

The output of the function is an object with properties corresponding to the following thermodynamic quantities
(all but the last three are from [lbftd](https://github.com/jmichaelb/LocalBasisFunction/tree/master/Python/lbftd)):

| Quantity  (PT only)      |  Symbol in SeaFreeze  |  Unit (SI)  |
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

| Quantity  (requires PTM)      |  Symbol in SeaFreeze  |  Unit (SI)  |
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

 **NaN values returned when out of parameterization boundaries.**

### Example

```python
import numpy as np
from seafreeze import seafreeze as sf

# list supported phases
sf.phases.keys()

# evaluate thermodynamics for ice VI at 900 MPa and 255 K
PT = np.empty((1,), dtype='object')
PT[0] = (900, 255)
out = sf.seafreeze(PT, 'VI')
# view a couple of the calculated thermodynamic quantities at this P and T
out.rho     # density
out.Vp      # compressional wave velocity

# evaluate thermodynamics for water at three separate PT conditions
PT = np.empty((3,), dtype='object')
PT[0] = (441.0858, 313.95)
PT[1] = (478.7415, 313.96)
PT[2] = (444.8285, 313.78)
out = sf.seafreeze(PT, 'water1')
# values for output fields correspond positionally to (P,T) tuples 
out.H       # enthalpy

# evaluate ice V thermodynamics at pressures 400-500 MPa and temperatures 240-250 K
P = np.arange(400, 501, 2)
T = np.arange(240, 250.1, 0.5)
PT = np.array([P, T], dtype='object')
out = sf.seafreeze(PT, 'V')
# rows in output correspond to pressures; columns to temperatures
out.A       # Helmholtz energy
out.shear   # shear modulus
```


## `seafreeze.whichphase`: determining the stable phase of water

### Usage
Seafreeze also includes a function to determine which of the *supported* phases is stable
under the given pressure and temperature conditions. 
The function `whichphase` has a single parameter, `PT`, 
which requires the same format as in the `seafreeze` function.

The output of the function is a [Numpy array](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html)
with an integer indicating the phase number corresponding to the `PT` input.  The phase number 0 means 
liquid water, phase number 1 means ice Ih, phase number 3 means ice III, etc.  Points outside the range
of all phases will return `numpy.nan`.
- for a list of scattered (P,T) conditions, each value corresponds to the same index in the input
- for a grid of PT conditions, each row corresponds to a pressure and each column to a temperature from the input.

`seafreeze.phasenum2phase` can be used to map output phase numbers to a phase.  
Each item in this dictionary has the phase number as its key and the phase as the value. 

### Example

```python
import numpy as np
from seafreeze.seafreeze import seafreeze as sf

# determine the phase of water at 900 MPa and 255 K
PT = np.empty((1,), dtype=object)
PT[0] = (900, 255)
out = sf.whichphase(PT)
# map to a phase using phasenum2phase
sf.phasenum2phase(out[0])


# determine phase for three separate (P,T) conditions
PT = np.empty((3,), dtype=object)
PT[0] = (100, 200)
PT[1] = (400, 250)
PT[2] = (1000, 300)
out = sf.whichphase(PT)
# show phase for each (P,T)
[(PT, sf.phasenum2phase(pn)) for (PT, pn) in zip(PT, out)]

# find the likely phases at pressures 0-5 MPa and temperatures 240-300 K
P = np.arange(0, 5, 0.1)
T = np.arange(240, 300)
PT = np.array([P, T], dtype=object)
out = sf.whichphase(PT)
```

## Important remarks 
### Water representation
The ices Gibbs parameterizations are optimized to be used with 'water1' Gibbs LBF from Bollengier et al. (2019), specially for phase equilibrium calculation. Using other water parameterization wil lead to incorrect melting curves. 'water2' (Brown 2018) and 'water_IAPWS95' (IAPWS95) parametrization are provided for HP extention (up to 100 GPa) and comparison only. The authors recommend the use of 'water1' (Bollengier et al. 2019) for any application in the 200-355 K range and up to 2300 MPa.

### Range of validity
SeaFreeze stability prediction is currently considered valid down to 130K, which correspond to the ice VI - ice XV transition. The ice Ih - II transition is potentially valid down to 73.4 K (ice Ih - ice XI transition). The ice VII and ice X representation extend to 1TPa (1e6 MPa) and 2000K.

## References
- [Bollengier, Brown and Shaw (2019) J. Chem. Phys. 151, 054501; doi: 10.1063/1.5097179](https://aip.scitation.org/doi/abs/10.1063/1.5097179)
- [Brown (2018) Fluid Phase Equilibria 463, pp. 18-31](https://www.sciencedirect.com/science/article/pii/S0378381218300530)
- [Feistel and Wagner (2006), J. Phys. Chem. Ref. Data 35, pp. 1021-1047](https://aip.scitation.org/doi/abs/10.1063/1.2183324)
- [Journaux et al., (2019), in review in JGR: Planets (available on ArXiv)](https://arxiv.org/abs/1907.09598)
- [Wagner and Pruss (2002), J. Phys. Chem. Ref. Data 31, pp. 387-535](https://aip.scitation.org/doi/abs/10.1063/1.1461829)
- [French and Redmer (2015), Physical Review B 91, 014308](http://link.aps.org/doi/10.1103/PhysRevB.91.014308)

## Authors

* **Baptiste Journaux** - *University of Washington, Earth and Space Sciences Department, Seattle, USA* 
* **J. Michael Brown** - *University of Washington, Earth and Space Sciences Department, Seattle, USA* 
* **Penny Espinoza** - *University of Washington, Earth and Space Sciences Department, Seattle, USA* 
* **Erica Clinton** - *University of Washington, Earth and Space Sciences Department, Seattle, USA* 
* **Tyler Gordon** - *University of Washington, Department of Astronomy, Seattle, USA*
* **Ula Jones** - *University of Washington, Earth and Space Sciences Department, Seattle, USA*

## Change log

### Changes since 0.9.0
- `1.0`: added NaCl aqueous solution EOS and concentration dependent thermodynamic variables.
- `0.9.4`: Adjusted python readme syntax and package authorship info
- `0.9.3`: add ice VII and ice X from French and Redmer (2015). LocalBasisFunction spline interpretation software integrated into SeaFreeze Python package. Adjusted packaging to work better with pip
- `0.9.2.post2`: `whichphase` returns `numpy.nan` if PT is outside the regime of all phases
- `0.9.2`: add ice II to the representation.
- `0.9.1`: add `whichphase` function

### Changes from 0.8
- rename function get_phase_thermodynamics to seafreeze
- reverse order of PT and phase in function signature
- remove a layer of nesting (`seafreeze.seafreeze` rather than `seafreeze.seafreeze.seafreeze`)


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

Illustration montage uses pictures from NASA Galileo and Cassini spacecrafts (from top to bottom: Enceladus, Europa and Ganymede). Terrestrial sea ice picture use with the authorization of the author [Rowan Romeyn](https://arcex.no/meet-rowan-romeyn-a-new-arcex-phd-student/).
