# SeaFreeze

V0.9.2

The SeaFreeze package allows computation of the thermodynamic and elastic properties of water and ice polymorphs Ih, III, V and VI in the 0-2300 MPa and 220-500 K range. It is based on the evaluation of Local Basis Functions for each phase. The formalism is described in more details in Brown (2018), Journaux et al. (2019), and in the liquid water Gibbs parameterization by Bollengier, Brown, and Shaw (2019).


## Installation
This package will install [uw-highP-geophysics-tools](https://pypi.org/project/uw-highP-geophysics-tools/) and its dependencies.

Run the following command to install

`pip3 install SeaFreeze`

To upgrade to the latest version, use

`pip3 install --upgrade SeaFreeze`


## `seafreeze.seafreeze`: calculating thermodynamic and elastic properties of a phase of water

### Usage
The main function of SeaFreeze is `seafreeze.seafreeze`, which has the following parameters:
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
  - 'water1' -  extends to 500 K and 2300 MPa; from Bollengier et al. 2019
  - 'water2' -  extends to 100 GPa; from Brown 2018
  - 'water_IAPWS95' - LBF representation of IAPWS 95; from Wagner and Pruß, 2002


The output of the function is an object with properties corresponding to the following thermodynamic quantities
(all but the last three are from [lbftd](https://github.com/jmichaelb/LocalBasisFunction/tree/master/Python/lbftd)):

| Quantity        |  Symbol in SeaFreeze  |  Unit (SI)  |
| --------------- |:---------------------:| :----------:|
| Gibbs Energy           | `G` | J/kg |
| Entropy                | `S` | J/K/kg |
| Internal Energy        | `U` | J/kg |
| Enthalpy               | `H` | J/kg |
| Helmholtz free energy  | `A` | J/kg |
| Density                |`rho`| kg/m<sup>3</sup> |
|Specific heat capacity at constant pressure|`Cp`| J/kg/K |
|Specific heat capacity at constant volume|`Cv`| J/kg/K |
| Isothermal bulk modulus      |`Kt`| MPa |
|Pressure derivative of the Isothermal bulk modulus|`Kp`| - |
| Isoentropic bulk modulus     |`Ks`| MPa |
| Thermal expansivity     |`alpha`| K<sup>-1</sup>  |
| Shear modulus     |`shear`| MPa |
| P wave velocity     |`Vp`| m/s |
| S wave velocity     |`Vs`| m/s |
| Bulk sound speed     |`vel`| m/s |

 **NaN values returned when out of parameterization boundaries.**

### Example

```python
import numpy as np
import seafreeze as sf

# list supported phases
sf.phases.keys()

# evaluate thermodynamics for ice VI at 900 MPa and 255 K
PT = np.empty((1,), np.object)
PT[0] = (900, 255)
out = sf.seafreeze(PT, 'VI')
# view a couple of the calculated thermodynamic quantities at this P and T
out.rho     # density
out.Vp      # compressional wave velocity

# evaluate thermodynamics for water at three separate PT conditions
PT = np.empty((3,), np.object)
PT[0] = (441.0858, 313.95)
PT[1] = (478.7415, 313.96)
PT[2] = (444.8285, 313.78)
out = sf.seafreeze(PT, 'water1')
# values for output fields correspond positionally to (P,T) tuples 
out.H       # enthalpy

# evaluate ice V thermodynamics at pressures 400-500 MPa and temperatures 240-250 K
P = np.arange(400, 501, 2)
T = np.arange(240, 250.1, 0.5)
PT = np.array([P, T])
out = sf.seafreeze(PT, 'V')
# rows in output correspond to pressures; columns to temperatures
out.A       # Helmholtz energy
out.shear   # shear modulus
```


## `seafreeze.whichphase`: determining the stable phase of water

### Usage
Seafreeze also includes a function to determine which of the *supported* phases is stable
under the given pressure and temperature conditions. 
The function `seafreeze.whichphase` has a single parameter, `PT`, 
which requires the same format as in the `seafreeze.seafreeze` function.

The output of the function is a [Numpy array](https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html)
with an integer indicating the phase number corresponding to the `PT` input.  The phase number 0 means 
liquid water, phase number 1 means ice Ih, phase number 3 means ice III, etc.
- for a list of scattered (P,T) conditions, each value corresponds to the same index in the input
- for a grid of PT conditions, each row corresponds to a pressure and each column to a temperature from the input.

`seafreeze.phasenum2phase` can be used to map output phase numbers to a phase.  
Each item in this dictionary has the phase number as its key and the phase as the value. 

### Example

```python
import numpy as np
import seafreeze as sf

# determine the phase of water at 900 MPa and 255 K
PT = np.empty((1,), np.object)
PT[0] = (900, 255)
out = sf.whichphase(PT)
# map to a phase using phasenum2phase
sf.phasenum2phase[out[0]]


# determine phase for three separate (P,T) conditions
PT = np.empty((3,), np.object)
PT[0] = (100, 200)
PT[1] = (400, 250)
PT[2] = (1000, 300)
out = sf.whichphase(PT)
# show phase for each (P,T)
[(PT, sf.phasenum2phase[pn]) for (PT, pn) in zip(PT, out)]

# find the likely phases at pressures 0-5 MPa and temperatures 240-300 K
P = np.arange(0, 5, 0.1)
T = np.arange(240, 300)
PT = np.array([P, T])
out = sf.whichphase(PT)
```

## Important remarks 
### Water representation
The ices Gibbs parameterizations are optimized to be used with 'water1' Gibbs LBF from Bollengier et al. (2019), specially for phase equilibrium calculation. Using other water parameterization wil lead to incorrect melting curves. 'water2' (Brown 2018) and 'water_IAPWS95' (IAPWS95) parametrization are provided for HP extention (up to 100 GPa) and comparison only. The authors recommend the use of 'water1' (Bollengier et al. 2019) for any application in the 200-355 K range and up to 2300 MPa.

### Range of validity
SeaFreeze stability prediction is currently considered valid down to 130K, which correspond to the ice VI - ice XV transition. The ice Ih - II transition is potentially valid down to 73.4 K (ice Ih - ice XI transition).

## References
- [Bollengier, Brown and Shaw (2019) J. Chem. Phys. 151, 054501; doi: 10.1063/1.5097179](https://aip.scitation.org/doi/abs/10.1063/1.5097179)
- [Brown (2018) Fluid Phase Equilibria 463, pp. 18-31](https://www.sciencedirect.com/science/article/pii/S0378381218300530)
- [Feistel and Wagner (2006), J. Phys. Chem. Ref. Data 35, pp. 1021-1047](https://aip.scitation.org/doi/abs/10.1063/1.2183324)
- [Journaux et al., (2019), in review in JGR: Planets (available on ArXiv)](https://arxiv.org/abs/1907.09598)
- [Wagner and Pruss (2002), J. Phys. Chem. Ref. Data 31, pp. 387-535](https://aip.scitation.org/doi/abs/10.1063/1.1461829)

## Author

* **Penny Espinoza** - *University of Washington, Earth and Space Sciences Department, Seattle, USA* 
* **Baptiste Journaux** - *University of Washington, Earth and Space Sciences Department, Seattle, USA* 
* **J. Michael Brown** - *University of Washington, Earth and Space Sciences Department, Seattle, USA* 

## Change log

### Changes since 0.9.0
- `0.9.1`: add `whichphase` function

### Changes from 0.8
- rename function get_phase_thermodynamics to seafreeze
- reverse order of PT and phase in function signature
- remove a layer of nesting (`seafreeze.seafreeze` rather than `seafreeze.seafreeze.seafreeze`)



