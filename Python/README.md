# SeaFreeze

V0.9.1

The SeaFreeze package allows computation of the thermodynamic and elastic properties of water and ice polymorphs Ih, III, V and VI in the 0-2300 MPa and 220-500 K range. It is based on the evaluation of Local Basis Functions for each phase. The formalism is described in more details in Brown (2018), Journaux et al. (2019), and in the liquid water Gibbs parametrization by Bollengier, Brown, and Shaw (2019).


## Installation
This package will install [uw-highP-geophysics-tools](https://pypi.org/project/uw-highP-geophysics-tools/) and its dependencies.

Run the following command to install

`pip3 install SeaFreeze`


## Usage

### `seafreeze.seafreeze`
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

### `seafreeze.whichphase`
Seafreeze also includes a function to determine which of the supported phases is most likely to exist
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





## Examples

### Calculating thermodynamics for a phase

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

### Determining the probable phase 

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

## Change log

### Changes from 0.8
- rename function get_phase_thermodynamics to seafreeze
- reverse order of PT and phase in function signature
- remove a layer of nesting (`seafreeze.seafreeze` rather than `seafreeze.seafreeze.seafreeze`)

### Changes since 0.9.0
- `0.9.1`: add `whichphase` function



