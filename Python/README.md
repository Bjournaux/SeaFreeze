# SeaFreeze

V0.8.0 (Python version)

The SeaFreeze package allows computation of the thermodynamic and elastic properties of water and ice polymorphs Ih, III, V and VI in the 0-2300 MPa and 220-500 K range. It is based on the evaluation of Local Basis Functions for each phase. The formalism is described in more details in Brown (2018), Journaux et al. (2019), and in the liquid water Gibbs parametrization by Bollengier, Brown, and Shaw (2019).


## Installation
This package will install [uw-highP-geophysics-tools](https://github.com/jmichaelb/LocalBasisFunction) and its dependencies.

Run the following command in Terminal on OSX or in Git Bash on Windows.

`pip3 install SeaFreeze`

Note that scripts can be run in the regular command prompt on Windows; it is only the install that needs to be done in Git Bash.

## Usage
The main function of SeaFreeze is `seafreeze.seafreeze.get_phase_thermodynamics`, which has the following parameters:
- `phase`: indicates the phase of H₂O.  Supported phases are
  - 'Ih' - from  Feistel and Wagner, 2006
  - 'III' - from Journaux et al, 2019
  - 'V' - from Journaux et al, 2019
  - 'VI' - from Journaux et al, 2019
  - 'water1' -  extends to 500 K and 2300 MPa; from Bollengier et al 2019
  - 'water2' -  extends to 100 GPa; from Brown 2018
  - 'water_IAPWS95' - LBF representation of IAPWS 95; from Wagner and Pruß, 2002
- `PT`: the pressure (MPa) and temperature (K) conditions at which the thermodynamic quantities should be
  calculated -- note that these are required units, as conversions are built into several calculations
  This parameter can have one of the following formats:
  - a 1-dimensional numpy array of tuples with one or more scattered (P,T) tuples 
  - a numpy array with 2 nested numpy arrays, the first with pressures and the second
    with temperatures -- each inner array must be sorted from low to high values
    a grid will be constructed from the P and T arrays such that each row of the output
    will correspond to a pressure and each column to a temperature 
- `path`: optional parameter with the location of the spline file to be evaluated.  The
    default value assumes the file is distributed with this project.

The output of the function is an object with properties corresponding to thermodynamic quantities, most
of which are listed in the [lbftd](https://github.com/jmichaelb/LocalBasisFunction/tree/master/Python/lbftd) README.
Also included are
- `Vp`: compressional wave velocity, in m/s
- `Vs`: shear wave velocity, in m/s
- `shear`: shear modulus, in MPa

```python
import numpy as np
import seafreeze.seafreeze as sf

# list supported phases
sf.phases.keys()

# evaluate thermodynamics for ice VI at 900 MPa and 255 K
PT = np.empty((1,), np.object)
PT[0] = (900, 255)
out = sf.get_phase_thermodynamics('VI', PT)
# show density and compressional wave velocity 
out.rho
out.Vp

# evaluate ice V thermodynamics for ice IV at three separate PT conditions
PT = np.empty((3,), np.object)
PT[0] = (441.0858, 313.95)
PT[1] = (478.7415, 313.96)
PT[2] = (444.8285, 313.78)
out = sf.get_phase_thermodynamics('V', PT)

# evaluate water's thermodynamics for a grid of pressures and temperatures
P = np.arange(0.1, 1000.2, 10)
T = np.arange(240, 501, 2)
PT = np.array([P, T])
out = sf.get_phase_thermodynamics('water1', PT)
```

