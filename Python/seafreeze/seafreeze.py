from collections import namedtuple
import os.path as op
import numpy as np
from mlbspline import load
from lbftd import evalGibbs as eg

defpath = op.join(op.dirname(op.abspath(__file__)), 'SeaFreeze_Gibbs.mat')

def get_phase_thermodynamics(phase, PT=None, path=defpath):
    """ Calculates thermodynamic quantities for H2O water or ice polymorphs Ih, III, V, and VI for all phases
        (see lbftd documentation for full list)
        for solid phases only:
            - Vp (compressional wave velocity, in m/s)
            - Vs (shear wave velocity, in m/s)
            - shear (shear modulus, in MPa)
    Requires the SeaFreeze_Gibbs.mat library containing the Gibbs LBF parametrization (installed with this module).

    NOTE:  The authors recommend the use of 'water1' for any application in the 200-355 K range and up to 2300 MPa.
    The ice Gibbs parametrizations are optimized for use with the 'water1' phase for phase equilibrium calculations.
    Using other water parametrizations will lead to incorrect melting curves -- 'water2' and 'water_IAPWS95'
    parametrizations are provided for HP extension up to 100 GPa and comparison only.

    :param phase:   one of the keys of the phases dict, indicating the phase of H2O to be evaluated
    :param PT:      the pressure (MPa) and temperature (K) conditions at which the thermodynamic quantities should be
                    calculated -- note that these are required units, as conversions are built into several calculations
                    This parameter can have one of the following formats:
                        - a 1-dimensional numpy array of tuples with one or more scattered (P,T) tuples, e.g.
                                PT = np.empty((3,), np.object)
                                PT[0] = (441.0858, 313.95)
                                PT[1] = (478.7415, 313.96)
                                PT[2] = (444.8285, 313.78)
                        - a numpy array with 2 nested numpy arrays, the first with pressures and the second
                          with temperatures -- each inner array must be sorted from low to high values
                          a grid will be constructed from the P and T arrays such that each row of the output
                          will correspond to a pressure and each column to a temperature, e.g.
                                P = np.arange(0.1, 1000.2, 10)
                                T = np.arange(240, 501, 2)
                                PT = np.array([P, T])
    :param path:    an optional path to the SeaFreeze_Gibbs.mat file
                    default path assumes the spline distributed along with the project
    :return:        object containing the calculated thermodynamic quantities (as named properties), as well as
                    a PTM property (a copy of PT)
    """
    try:
        phasedesc = phases[phase]
    except KeyError:
        raise ValueError('The specified phase is not recognized.  Supported phases are ' +
                         ', '.join(phases.keys())+'.')
    sp = load.loadSpline(path, phasedesc.sp_name)
    # calc density and isentropic bulk modulus
    isscatter = _is_scatter(PT)
    tdvs = _get_tdvs(sp, PT, isscatter)
    if phasedesc.shear_mod_parms:
        smg = _get_shear_mod_GPa(phasedesc.shear_mod_parms, tdvs.rho, _get_T(PT, isscatter))
        tdvs.shear = 1e3 * smg  # convert to MPa for consistency with other measures
        tdvs.Vp = _get_Vp(smg, tdvs.rho, tdvs.Ks)
        tdvs.Vs = _get_Vs(smg, tdvs.rho)
    return tdvs


def _get_tdvs(sp, PT, is_scatter):
    """ peeks into PT to see if the PT data is for grid or scatter and calls the appropriate evalGibbs function

    :param sp:  the Gibbs LBF
    :param PT:  the PT data
    :return:    tdv object
    """
    fn = eg.evalSolutionGibbsScatter if is_scatter else eg.evalSolutionGibbsGrid
    return fn(sp, PT, failOnExtrapolate=False)


def _get_shear_mod_GPa(sm, rho, T):
    return None if sm is None else sm[0] + sm[1]*(rho - sm[4]) + sm[2]*(rho-sm[4])**2 + sm[3]*(T-sm[5])


def _get_Vp(smg, rho, Ks):
    return 1e3 * np.sqrt((Ks/1e3 + 4/3*smg)/rho/1e-3)


def _get_Vs(smg, rho):
    return 1e3 * np.sqrt(smg/rho/1e-3)


def _is_scatter(PT):
    return isinstance(PT[0], tuple) or (PT.shape == (1,2) and np.isscalar(PT[0]) and np.isscalar(PT[1]))


def _get_T(PT, is_scatter):
    return np.array([T for P,T in PT]) if is_scatter else PT[1]


#########################################
## Constants
#########################################
PhaseDesc = namedtuple('PhaseDesc', 'sp_name shear_mod_parms')
defpath = op.join(op.dirname(op.abspath(__file__)), 'SeaFreeze_Gibbs.mat')
phases = {"Ih": PhaseDesc("G_iceIh", [3.04, -0.00462, 0, -0.00607, 1000, 273.15]),  # Feistel and Wagner, 2006
          "III": PhaseDesc("G_iceIII", [2.57, 0.0175, 0, -0.014, 1100, 273]),       # Journaux et al, 2019
          "V": PhaseDesc("G_iceV", [2.57, 0.0175, 0, -0.014, 1100, 273]),           # Journaux et al, 2019
          "VI": PhaseDesc("G_iceVI", [2.57, 0.0175, 0, -0.014, 1100, 273]),         # Journaux et al, 2019
          "water1": PhaseDesc("G_H2O_2GPa_500K", None),     # extends to 500 K and 2300 MPa; Bollengier et al 2019
          "water2": PhaseDesc("G_H2O_100GPa_10000K", None), # extends to 100 GPa; Brown 2018
          "water_IAPWS95": PhaseDesc("G_H2O_IAPWS", None)   # LBF representation of IAPWS 95; Wagner and Pruß, 2002
          }
