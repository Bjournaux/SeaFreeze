from collections.__init__ import namedtuple
from inspect import getargvalues, currentframe
from pprint import pformat

import numpy as np
from numpy.lib.scimath import sqrt

from mlbspline.eval import evalMultivarSpline


'''
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Warning: units must be as specified in the README because some conversions are hardcoded.
With the exception of pressure, units are SI.  Pressure is in MPa rather than Pa.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'''


#########################################
## TDV eval functions
#########################################

def evalGibbsEnergy(gibbsSp, PTM):
    """
    :return: G
    """
    return evalMultivarSpline(gibbsSp, PTM)


def evalIsobaricSpecificHeat(gibbsSp, gPTM, derivs):
    """
    :return:        Cp
    """
    if gibbsSp['ndT']:
        T = np.exp(gPTM[iT])*gibbsSp['Tc']
    else:
        T = gPTM[iT]
    return -1 * derivs.d2T * T


def evalIsochoricSpecificHeat(tdv, gPTM, derivs):
    """

    :return:    Cv
    """
    #Cv= Cp + Tm.*dPT.^2./d2P;
    return tdv.Cp + gPTM[iT] * np.power(derivs.dPT, 2) / derivs.d2P


def evalEntropy(derivs):
    """
    :return:        S
    """
    return -1 * derivs.d1T


def evalSoundSpeed(derivs):
    """
    :return:        vel
    """
    # MPa-Pa units conversion cancels
    return np.real(sqrt(np.power(derivs.d1P, 2) / (np.power(derivs.dPT, 2) / derivs.d2T - derivs.d2P)))


def evalDensity(derivs):
    """
    :return:        rho
    """
    return 1e6 * np.power(derivs.d1P, -1)   # 1e6 for MPa to Pa


def evalIsothermalBulkModulus(derivs):
    """
    :return:        Kt
    """
    return -1 * derivs.d1P / derivs.d2P


def evalPDerivIsothermalBulkModulus(derivs):
    """
    :return:        Kp (K')
    """
    return derivs.d1P * np.power(derivs.d2P, -2) * derivs.d3P - 1


def evalIsentropicBulkModulus(tdv):
    """
    :return:        Ks
    """
    return tdv.rho * np.power(tdv.vel, 2) / 1e6 #  1e6 for MPa conversion


def evalThermalExpansivity(derivs, tdv):
    """
    :return:        alpha
    """
    return 1e-6 * derivs.dPT * tdv.rho  #  1e-6 for MPa conversion


def evalInternalEnergy(gPTM, tdv):
    """
    :param gPTM:    gridded dimensions over which spline is being evaluated
    :return:        U
    """
    return tdv.G - 1e6 * gPTM[iP] / tdv.rho + gPTM[iT] * tdv.S  #  1e6 for MPa conversion


def evalSoluteChemicalPotential(MWu, f, derivs, tdv):
    """
    :return:        mus
    """
    return MWu * tdv.G + f * derivs.d1M


def evalSolventChemicalPotential(MWv, f, gPTM, derivs, tdv):
    """
    :return:        muw
    """
    return (tdv.G * MWv) - (MWv * f * gPTM[iM] * derivs.d1M)


def evalHelmholtzEnergy(gPTM, tdv):
    """
    :return:        A
    """
    return tdv.U - gPTM[iT] * tdv.S


def evalEnthalpy(gPTM, tdv):
    """

    :return:        H
    """
    return tdv.G + gPTM[iT] * tdv.S


def evalPartialMolarVolume(MWu, f, derivs):
    """ Slope at a point of the V v. M graph (volume vs molality)
    :return:        Vm
    """
    return (MWu * derivs.d1P) + (f * derivs.dPM)


def evalPartialMolarHeatCapacity(MWu, f, tdv, derivs, gPTM):
    """
    :return:    Cpm
    """
    return MWu * tdv.Cp - f * derivs.d2T1M * gPTM[iT]


def evalVolume(tdv):
    """
    :return:    V
    """
    return np.power(tdv.rho, -1)


def evalApparentSpecificHeat(f, gPTM, tdv):
    """
    :return:    Cpa
    """
    return (tdv.Cp * f - _get0MTdv(tdv, 'Cp')) / _getDividableBy(gPTM[iM])


def evalApparentVolume(f, gPTM, tdv):
    """ slope of a chord between pure solvent and a concentration on a V v. M graph
    :return:    Va
    """
    return 1e6 * (tdv.V * f - _get0MTdv(tdv, 'V')) / _getDividableBy(gPTM[iM]) # 1e6 for MPa conversion


def _get0MTdv(tdv, prop):
    slc = [slice(None)] * (len(tdv.PTM))
    slc[iM] = slice(0, 1)
    return getattr(tdv, prop)[tuple(slc)]


def _getDividableBy(inp):
    eps = np.finfo(inp.dtype).eps
    return np.where(inp != 0, inp, eps)


#########################################
## Set up derivatives
#########################################

def _getSupportedDerivatives():
    """ Define a name for a derivative of Gibbs energy
    and the derivative order with respect to P, T, and M required to calculate it
    (if not provided, all default to defDer)
    These can then be used in a TDVSpec

    :return:    an immutable iterable containing a list of supported derivatives
    """
    return tuple([
        derivSpec('d1P', wrtP=1, wrtT=defDer, wrtM=defDer),
        derivSpec('d1T', wrtP=defDer, wrtT=1, wrtM=defDer),
        derivSpec('d1M', wrtP=defDer, wrtT=defDer, wrtM=1),
        derivSpec('dPT', wrtP=1, wrtT=1, wrtM=defDer),
        derivSpec('dPM', wrtP=1, wrtT=defDer, wrtM=1),
        derivSpec('d2P', wrtP=2, wrtT=defDer, wrtM=defDer),
        derivSpec('d2T', wrtP=defDer, wrtT=2, wrtM=defDer),
        derivSpec('d2T1M', wrtP=defDer, wrtT=2, wrtM=1),
        derivSpec('d3P', wrtP=3,wrtT=defDer, wrtM=defDer)
    ])


#########################################
## Set up thermodynamic variables (TDVs)
#########################################

def _getTDVSpec(name, calcFn, reqM=False, reqMWv=False, parmMWv='MWv', reqMWu=False, parmMWu='MWu',
                reqGrid=False, parmgrid='gPTM', reqF=False, parmf='f',
                reqDerivs=[], parmderivs='derivs', reqTDV=[], parmtdv='tdv', reqSpline=False,
                parmspline='gibbsSp', reqPTM=False, parmptm='PTM', req0M=False):
    """ Builds a TDVSpec namedtuple indicating what is required to calculate this particular thermodynamic variable
    :param name:        the name / symbol of the tdv (e.g., G, rho, alpha, muw)
    :param calcFn:      the name of the function used to calculate the tdv
    :param reqM:        True if DIRECTLY calculating the tdv requires concentration.  False otherwise
    :param reqMWv:      True if DIRECTLY calculating the tdv requires solvent molecular weight.  False otherwise
    :param parmMWv:     the name of the parameter of calcFn used to pass in solvent molecular weight if reqMWv
    :param reqMWu:      True if DIRECTLY calculating the tdv requires solute molecular weight.  False otherwise
    :param parmMWu:     the name of the parameter of calcFn used to pass in solute molecular weight if reqMWu
    :param reqGrid:     True if DIRECTLY calculating the tdv requires you to grid the PT[M] values.  False otherwise
    :param parmgrid:    the name of the parameter of calcFn used to pass in the gridded input dimensions if reqGrid
    :param reqF:        True if DIRECTLY calculating the tdv requires a conversion factor that gives the volume
                        of solvent in a unit volume of solution (this is calculated by
                        evalGibbs._getVolSolventInVolSlnConversion)
    :param parmf:       the name of the parameter of calcFn used to pass in the vol solvent to vol solution conversion
                        factor
    :param reqDerivs:   A list of derivatives needed to DIRECTLY calculate the tdv
                        (e.g. for tdv Kp, this would be ['d1P','dP','d3P']
                        - see fn evalPDerivIsothermalBulkModulus)
                        See getDerivatives for a full list of derivatives that can be calculated
    :param parmderivs:  the name of the parameter of calcFn used to pass in the pre-calculated derivatives if reqDerivs
    :param reqTDV:      A list of other thermodynamic variables needed to DIRECTLY calculate the tdv
                        (e.g. for tdv U, this would be ['G','rho'] - see fn evalInternalEnergy)
                        Note: recursive dependencies are not supported
    :param parmtdv:     the name of the parameter of calcFn used to pass in the tdvout if reqTDV
    :param reqSpline:   If True, calcFn needs the spline definition
    :param parmspline:  the name of the parameter of calcFn used to pass in the spline def if reqSpline
    :param reqPTM:      If True, calcFn needs the original dimension input (parm PTM in evalSolutionGibbs*) to run
    :param parmptm:     the name of the parameter of calcFn used to pass in the original input if reqPTM
    :param req0M:       if True, calcFn needs the 0 concentration for calculating apparent values
    :return:            a namedtuple giving the spec for the tdv
    """
    reqTDV = set(reqTDV); reqDerivs = set(reqDerivs)    # make these sets to enforce uniqueness and immutability
    if name in reqTDV:
        raise ValueError('Recursive dependencies are not supported.  Amend the ' + name + ' TDVSpec accordingly')
    if reqDerivs.difference(derivativenames):
        raise ValueError('One or more derivatives are not supported. Amend the ' + name + ' TDVSpec accordingly.' +
                         'Supported derivatives are '+pformat(derivativenames))
    arginfo = getargvalues(currentframe())
    # build properties of the TDVSpec dynamically
    flds = arginfo.args
    if not reqMWv:      flds.remove('parmMWv')
    if not reqMWu:      flds.remove('parmMWu')
    if not reqGrid:     flds.remove('parmgrid')
    if not reqDerivs:   flds.remove('parmderivs')
    if not reqTDV:      flds.remove('parmtdv')
    if not reqSpline:   flds.remove('parmspline')
    if not reqPTM:      flds.remove('parmptm')
    vals = {f: v for f, v in arginfo.locals.items() if f in flds}

    tdvspec = namedtuple('TDVSpec', flds)
    return tdvspec(**vals)


def _getSupportedThermodynamicVariables():
    """ When adding new tdvs, you don't need to worry about adding them in any particular order -
        dependencies are handled elsewhere
        See also the comments for _getTDVSpec and evalSolutionGibbs*

    :return: immutable iterable with the the full set of specs for supported thermodynamic variables
    """
    out = tuple([
        _getTDVSpec('G', evalGibbsEnergy, reqSpline=True, reqPTM=True),
        _getTDVSpec('U', evalInternalEnergy, reqGrid=True, reqTDV=['G', 'rho', 'S']),
        _getTDVSpec('A', evalHelmholtzEnergy, reqGrid=True, reqTDV=['U', 'S']),
        _getTDVSpec('H', evalEnthalpy, reqGrid=True, reqTDV=['G', 'S']),
        _getTDVSpec('S', evalEntropy, reqDerivs=['d1T']),
        _getTDVSpec('rho', evalDensity, reqDerivs=['d1P']),
        _getTDVSpec('V', evalVolume, reqTDV=['rho']),
        _getTDVSpec('Cp', evalIsobaricSpecificHeat, reqSpline=True, reqGrid=True, reqDerivs=['d2T']),
        _getTDVSpec('Cv', evalIsochoricSpecificHeat, reqGrid=True, reqDerivs=['dPT', 'd2P'], reqTDV=['Cp']),
        _getTDVSpec('Kt', evalIsothermalBulkModulus, reqDerivs=['d1P', 'd2P']),
        _getTDVSpec('Kp', evalPDerivIsothermalBulkModulus, reqDerivs=['d1P', 'd2P', 'd3P']),
        _getTDVSpec('Ks', evalIsentropicBulkModulus, reqTDV=['rho', 'vel']),
        _getTDVSpec('vel', evalSoundSpeed, reqDerivs=['d1P', 'dPT', 'd2T', 'd2P']),
        _getTDVSpec('alpha', evalThermalExpansivity, reqDerivs=['dPT'], reqTDV=['rho']),
        _getTDVSpec('mus', evalSoluteChemicalPotential, reqMWu=True, reqF=True, reqDerivs=['d1M'], reqTDV=['G']),
        _getTDVSpec('muw', evalSolventChemicalPotential, reqM=True, reqMWv=True, reqGrid=True, reqF=True,
                    reqDerivs=['d1M'], reqTDV=['G']),
        _getTDVSpec('Vm', evalPartialMolarVolume, reqMWu=True, reqF=True, reqDerivs=['d1P', 'dPM']),
        _getTDVSpec('Cpm', evalPartialMolarHeatCapacity, reqMWu=True, reqGrid=True, reqF=True, reqDerivs=['d2T1M'],
                    reqTDV=['Cp']),
        _getTDVSpec('Cpa', evalApparentSpecificHeat, reqM=True, reqGrid=True, reqF=True, reqTDV=['Cp'], req0M=True),
        _getTDVSpec('Va', evalApparentVolume, reqM=True, reqGrid=True, reqF=True, reqTDV=['V'], req0M=True)
            ])
    # check that all reqTDVs are represented in the list
    outnames = frozenset([t.name for t in out])
    unrecognizedDependencies = [(tdv.name, dep) for tdv in out for dep in tdv.reqTDV if dep not in outnames]
    if unrecognizedDependencies:
        raise ValueError('One or more statevars depend on an unrecognized TDV: '+pformat(unrecognizedDependencies))
    return out, outnames


def _getPTOnlyTDVSpec():
    """ determines the limited TDVs that can be calculated for a 2-D (P and T only) spline
    """
    # concentration required if the flag says it does, if it uses F (which req concentration),
    # or if it uses a concentration derivative
    Mreq = [m for m in statevars if m.reqM or m.reqF or [d for d in m.reqDerivs if 'M' in d]]
    # return values should be immutable - use tuples as TDVSpec namedtuples are not hashable so set is not usable
    return tuple([m for m in statevars if m not in Mreq])


#########################################
## Interpret and expand a spec
#########################################

def _addTDVDependencies(origTDVs):
    """ recursively adds thermodynamic variables on which origTDVs are dependent
    This does not handle derivative dependencies - see getDerivatives

    :param origTDVs:    an iterable of TDVs to check for dependencies
                        elements can be either the string names of statevars or the TDV objects
    :return:            immutable iterable of TDV objects including origTDVs and any TDVs on which they are
                        directly or indirectly dependent
    """
    otdvnames = set(origTDVs) if isinstance(next(iter(origTDVs)), str) else {m.name for m in origTDVs}
    while True:
        # get TDVs dependent on TDVs in otdvnames
        reqot = {t for m in statevars if m.name in otdvnames for t in m.reqTDV}
        # if every tdv in reqot is already in otdvnames, stop
        if reqot.issubset(otdvnames):
            break
        else:
            # add the new list of tdvs and repeat the loop to add more nested dependencies
            otdvnames = otdvnames.union(reqot)
    return tuple([m for m in statevars if m.name in otdvnames])


def expandTDVSpec(tdvSpec, dimCt):
    """ add dependencies for requested thermodynamic variables or sets default tdvSpec if none provided
    derivatives are handled separately - see getDerivatives

    :param tdvSpec: an iterable with the names of thermodynamic variables to be evaluated
    :param dimCt:   the number of dimensions - used if setting the default tdvSpec
    :return:        an immutable iterable of TDV namedtuples that includes those specified by tdvSpec
                    and all their dependencies
    """
    # if no spec provided, use the default spec (based on # dims of the spline)
    if len(tdvSpec) == 0:
        if dimCt == 3:
            tdvSpec = list(statevarnames)   # make mutable
        else:
            tdvSpec = [t.name for t in tdvsPTOnly]
    # TODO: try to get rid of this (horrible) next line
    tdvSpec = tdvSpec if not isinstance(tdvSpec, str) else (tdvSpec,)
    # check for completely unsupported statevars so no one has to evaluate twice because of a typo
    unsupported = set(tdvSpec) - statevarnames
    if unsupported:
        raise ValueError('One or more unsupported statevars have been requested: ' + pformat(unsupported))

    tdvSpec = _addTDVDependencies(tdvSpec)  # add thermodynamic variables on which requested ones depend
    return tdvSpec


#########################################
## Constants
#########################################
iP = 0; iT = 1; iM = 2      # dimension indices
defDer = 0                  # default derivative

derivSpec = namedtuple('DerivativeSpec', ['name', 'wrtP', 'wrtT', 'wrtM'])
derivatives = _getSupportedDerivatives()
derivativenames = {d.name for d in derivatives}

statevars, statevarnames = _getSupportedThermodynamicVariables()
# [(sv.name, sv.calcFn.__name__[4:]) for sv in statevars]        # lists TDV symbols and names

tdvsPTOnly = _getPTOnlyTDVSpec()
