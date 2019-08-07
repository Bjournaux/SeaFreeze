# SeaFreeze

V0.9

The SeaFreeze package allows computation of the thermodynamic and elastic properties of water and ice polymorphs Ih, III, V and VI in the 0-2300 MPa and 220-500 K range. It is based on the evaluation of Local Basis Functions for each phase. The formalism is described in more details in Brown (2018), Journaux et al. (2019), and in the liquid water Gibbs parametrization by Bollengier, Brown, and Shaw (2019).


## Installation
This package will install [uw-highP-geophysics-tools](https://github.com/jmichaelb/LocalBasisFunction) and its dependencies.

Run the following command in Terminal on OSX or in Git Bash on Windows.

`pip3 install SeaFreeze`

Note that scripts can be run in the regular command prompt on Windows; it is only the install that needs to be done in Git Bash.

## Usage
The main function of SeaFreeze is `seafreeze.seafreeze.get_phase_thermodynamics`, which has the following parameters:
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
  - 'III' - from Journaux et al, 2019
  - 'V' - from Journaux et al, 2019
  - 'VI' - from Journaux et al, 2019
  - 'water1' -  extends to 500 K and 2300 MPa; from Bollengier et al 2019
  - 'water2' -  extends to 100 GPa; from Brown 2018
  - 'water_IAPWS95' - LBF representation of IAPWS 95; from Wagner and Pruß, 2002


The output of the function is an object with properties corresponding to thermodynamic quantities, most
of which are listed in the [lbftd](https://github.com/jmichaelb/LocalBasisFunction/tree/master/Python/lbftd) README.
Also included are
- `Vp`: compressional wave velocity, in m/s
- `Vs`: shear wave velocity, in m/s
- `shear`: shear modulus, in MPa

## Examples

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

## Changes from 0.8
- remove a layer of nesting
- rename function get_phase_thermodynamics to seafreeze
- reverse order of PT and phase in function signature

