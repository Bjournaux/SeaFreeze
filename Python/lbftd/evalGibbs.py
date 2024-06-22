from datetime import datetime
from pprint import pformat
from time import time
import logging
from collections.abc import Iterable
import numpy as np
from psutil import virtual_memory
from mlbspline.eval import evalMultivarSpline, _isExtrapolation
from lbftd import statevars
from lbftd.statevars import iT, iP, iM, _getSupportedThermodynamicVariables, eps
from inspect import getmembers

log = logging.getLogger('lbftd')
stream = logging.StreamHandler()
stream.setFormatter(logging.Formatter('[LBFTD %(levelname)s] %(message)s'))


def evalSolutionGibbsGrid(gibbsSp, PTM, *tdvSpec, verbose=False, failOnExtrapolate=False, allowExtrapolations=False):
    """ Calculates thermodynamic variables for solutions based on a spline giving Gibbs energy
    This currently only supports single-solute solutions.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Warning: units must be as specified in the README because some conversions are hardcoded into this function.
    With the exception of pressure, units are SI.  Pressure is in MPa rather than Pa.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    :param gibbsSp: A B-spline  (in format given by loadSpline.loadMatSpline) for giving Gibbs energy (J/mg)
                    with dimensions pressure (MPa), temperature (K), and (optionally) molality (mol/kg), in the
                    order defined by statevars.iP, iT, and iM.
                    If molality is not provided, this function assumes that it is calculating
                    thermodynamic properties for a pure substance.
                    To calculate some tdvs, this spline must also include the following values:
                        - MWv (molecular weight of solvent in kg/mol), defaults to molecular weight of water
                        - MWu (molecular weight of solute in kg/mol)
                        - nu (number of ions in solution, dimensionless)
    :param PTM:     a Numpy ndarray of ndarrays with the conditions at which gibbsSp should be evaluated
                    the number of dimensions must be same as in the spline (PTM.size == gibbsSp['number'].size)
                    each of the inner ndarrays represents one of pressure (P), temperature (T), or molality (M)
                    and must be in the same order and units described in the README.
                    Additionally, each dimension must be sorted from low to high values.
    :param tdvSpec: iterable indicating the thermodynamic variables to be calculated
                    elements can be either strings showing the names (full list at statevars.statevarnames)
                    or TDV objects from statevars.statevars.
                    If not provided, this function will calculate the variables in statevars.tdvsPTOnly for a PT spline,
                    and those in statevars.statevars for a PTM spline.
                    Any other args provided will result in an error.
                    See the README for units.
    :param verbose: bool indicating whether to print status updates, warnings, etc. Default to False.
    :param failOnExtrapolate:   boolean indicating whether to error out if PTM includes values that fall outside the
                    knot range of gibbsSp.  If False (default), proceeds with the calculation.
                    See also allowExtrapolations parameter.
    :param allowExtrapolations: boolean indicating whether to return calculated output for values of x that fall outside
                    the knot range of the spline. False (default) returns numpy.nan for such values of x.
                    See also failOnExtrapolate.
    :return:        a named tuple with the requested thermodynamic variables as named properties
                    matching the statevars requested in the *tdvSpec parameter of this function
                    the output will also include P, T, and M (if provided) properties
    """
    statevarnames = _getSupportedThermodynamicVariables()
    [dimCt, tdvSpec, addedTdvs] = _parseInput(gibbsSp, *tdvSpec)
    _checkInputs(gibbsSp, dimCt, tdvSpec, PTM, failOnExtrapolate)
    tdvout = createThermodynamicStatesObj(tdvSpec, PTM, initializetdvs=True)  # this has the original PTM
    if _needsCutoff(dimCt, tdvSpec):
        PTM = _getPTMWithCutoff(PTM, gibbsSp['cutoff'])
    tdvout = _evalInternal(gibbsSp, tdvSpec, PTM, tdvout, dimCt, allowExtrapolations, verbose)
    return tdvout


def evalSolutionGibbsScatter(gibbsSp, PTM, *tdvSpec, failOnExtrapolate=False,
                             verbose=False, allowExtrapolations=False):
    """ Calculates thermodynamic variables for solutions based on a spline giving Gibbs energy
    This currently only supports single-solute solutions.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Warning: units must be as specified in the README because some conversions are hardcoded into this function.
    With the exception of pressure, units are SI.  Pressure is in MPa rather than Pa.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    :param gibbsSp: A B-spline  (in format given by loadSpline.loadMatSpline) for giving Gibbs energy (J/mg)
                    with dimensions pressure (MPa), temperature (K), and (optionally) molality (mol/kg), in the
                    order defined by statevars.iP, iT, and iM.
                    If molality is not provided, this function assumes that it is calculating
                    thermodynamic properties for a pure substance.
                    To calculate some tdvs, this spline must include the following values:
                        - MWv (molecular weight of solvent in kg/mol), defaults to molecular weight of water
                        - MWu (molecular weight of solute in kg/mol)
                        - nu (number of ions in solution, dimensionless)
    :param PTM:     a Numpy ndarray of tuples with the points at which gibbsSp should be evaluated
                    the number of elements in each tuple must be same as in the spline (PTM[i].size == gibbsSp['number'].size)
                    each tuples gives the pressure (P), temperature (T), and (optionally) molality (M)
                    and must be in the same order and units described in the notes for the gibbsSp parameter
                    No sorting is required.
    :param tdvSpec: iterable indicating the thermodynamic variables to be calculated
                    elements can be either strings showing the names (full list at statevars.statevarnames)
                    or TDV objects from statevars.statevars.
                    If not provided, this function will calculate the variables in statevars.tdvsPTOnly for a PT spline,
                    and those in statevars.statevars for a PTM spline.
                    Any other args provided will result in an error.
                    See the README for units.
    :param failOnExtrapolate:   boolean indicating whether to error out if PTM includes values that fall outside the
                    knot range of gibbsSp.  If False (default), proceeds with the calculation.
                    See also allowExtrapolations parameter.
    :param allowExtrapolations: boolean indicating whether to return calculated output for values of x that fall outside
                    the knot range of the spline. False (default) returns numpy.nan for such values of x.
                    See also failOnExtrapolate.
    :param verbose: boolean indicating whether to print status updates, warnings, etc.
    :return:        a named tuple with the requested thermodynamic variables as named properties
                    matching the statevars requested in the *tdvSpec parameter of this function
                    Each statevar property will contain an ndarray of values corresponding to the points listed in PTM.
                    The output also contains a PTM property for convenience.
    """

    [dimCt, tdvSpec, addedTdvs] = _parseInput(gibbsSp, *tdvSpec)
    needsCutoff = _needsCutoff(dimCt, tdvSpec)
    fakePTM = _makeFakePTMGrid(dimCt, PTM)  # Fakes grid-style input to do all input checks before processing
    _checkInputs(gibbsSp, dimCt, tdvSpec, fakePTM, failOnExtrapolate)
    tdvout = createThermodynamicStatesObj(tdvSpec, PTM, initializetdvs=True)
    # Point by point
    for p in range(PTM.size):
        # this function will always replace 0 concentration with eps
        wPTM = _ptmTuple2NestedArrays(PTM[p], dimCt, needsCutoff)
        tempout = createThermodynamicStatesObj(tdvSpec, wPTM)
        if allowExtrapolations or not _isPointExtrapolation(dimCt, gibbsSp['knots'], wPTM):
            # Prepend a "cutoff" concentration if one is needed by any of the quantities being calculated
            # This is slow, but we expect scatter to only be called for a limited number of points
            if needsCutoff:
                wPTM = _getPTMWithCutoff(wPTM, gibbsSp['cutoff'])
            tempout = _evalInternal(gibbsSp, tdvSpec, wPTM, tempout, dimCt, allowExtrapolations, verbose)
            for t in tdvSpec:
                tolastpt = getattr(tdvout, t.name)
                tolastpt[p] = getattr(tempout, t.name).flat[0]
                setattr(tdvout, t.name, tolastpt)
    return tdvout


def _evalInternal(gibbsSp, tdvSpec, PTM, tdvout, dimCt, allowExtrapolations, verbose=False):
    derivs = getDerivatives(gibbsSp, PTM, dimCt, tdvSpec, allowExtrapolations, verbose)
    gPTM = _getGriddedPTM(tdvSpec, PTM, verbose) if any([tdv.reqGrid for tdv in tdvSpec]) else None
    f = _getVolSolventInVolSlnConversion(gibbsSp['MW'][1], PTM) if any([tdv.reqF for tdv in tdvSpec]) else None
    strip_cutoff = False

    # Calculate thermodynamic variables and store in appropriate fields in tdvout
    completedTDVs = set()  # list of completed tdvs
    while len(completedTDVs) < len(tdvSpec):
        # Get tdvs that either have empty reqTDV or all of those tdvs have already been calculated
        # but that are not in completedTDVs themselves (don't recalculate anything)
        tdvsToEval = tuple(t for t in tdvSpec if
                           t.name not in completedTDVs and (not t.reqTDV or not t.reqTDV.difference(completedTDVs)))
        for t in tdvsToEval:
            if t.reqCutoff:
                strip_cutoff = True
            args = _buildEvalArgs(t, derivs, PTM, gPTM, tdvout, gibbsSp, f, allowExtrapolations)
            start = time()
            setattr(tdvout, t.name, t.calcFn(**args))  # Calculate the value and set it in the output
            end = time()
            if verbose: _printTiming('tdv ' + t.name, start, end)
            completedTDVs.add(t.name)
    if strip_cutoff:
        _removeCutoff(tdvout)
    return tdvout



def _parseInput(gibbsSp, *tdvSpec):
    dimCt = gibbsSp['number'].size
    # Expand spec to add dependencies (or set to default spec if no spec given)
    origSpec = tdvSpec
    tdvSpec = statevars.expandTDVSpec(tdvSpec, dimCt)
    addedTDVs = [s.name for s in tdvSpec if s.name not in origSpec]
    if origSpec and addedTDVs:  # The original spec was not empty and more tdvs were added
        print('NOTE: The requested thermodynamic variables depend on the following variables, which will be ' +
              'included as properties of the output object: ' + pformat(addedTDVs))
    return dimCt, tdvSpec, addedTDVs


def _buildEvalArgs(tdv, derivs, PTM, gPTM, tdvout, gibbsSp, f, allowExtrapolations):
    # check inputs should have found any bad inputs before this is called
    args = dict()
    if tdv.reqDerivs: args[tdv.parmderivs] = derivs
    if tdv.reqGrid:   args[tdv.parmgrid] = gPTM
    if tdv.reqMWv:    args[tdv.parmMWv] = H20_MWv if 'MW' not in gibbsSp else gibbsSp['MW'][0]  # loadgibbs always wraps this in an array
    if tdv.reqMWu:    args[tdv.parmMWu] = gibbsSp['MW'][1] # currently only single solute solutions are supported
    if tdv.reqNu:     args[tdv.parmNu] = gibbsSp['nu']
    if tdv.reqTDV:    args[tdv.parmtdv] = tdvout
    if tdv.reqSpline: args[tdv.parmspline] = gibbsSp
    if tdv.reqPTM:    args[tdv.parmptm] = PTM
    if tdv.reqF:      args[tdv.parmf] = f
    if tdv.reqAllowExtrap: args[tdv.parmAllowExtrap] = allowExtrapolations
    return args


def _makeFakePTMGrid(dimCt, PTM):
    out = np.empty(dimCt, object)
    for d in np.arange(0, dimCt):
        u = {t[d] for t in PTM}  # Get unique values in each dimension
        out[d] = np.array(sorted(u))  # Sort the values in ascending order
    return out


def _checkInputs(gibbsSp, dimCt, tdvSpec, PTM, failOnExtrapolate):
    """ Checks error conditions before performing any calculations, throwing an error if anything doesn't match up
      - Ensures necessary data is available for requested statevars (check req parameters for statevars)
      - Throws a warning (or error if failOnExtrapolate=True) if spline is being evaluated on values
        outside the range of its knots
      - Warns the user if the requested output would exceed vmWarningFactor times the amount of virtual memory
      Note that the mlbspline.eval module performs additional checks (like dimension count mismatches
    """
    knotranges = [(k[0], k[-1]) for k in gibbsSp['knots']]
    # If a PTM spline, make sure the spline's concentration dimension starts with 0 if any tdv has req0M=True
    # Note that the evalGibbs functions elsewhere handle the case where PTM does not include 0 concentration.
    reqCutoff = [t for t in tdvSpec if t.reqCutoff]
    cutoff = None if 'cutoff' not in gibbsSp else gibbsSp['cutoff']
    if dimCt == 3 and reqCutoff and knotranges[iM][0] > cutoff: # 1st knot in concentration will no longer be 0
        raise ValueError('You cannot calculate ' + pformat([t.name for t in reqCutoff]) + ' with a spline that does ' +
                         'not include a concentration of ' + str(cutoff) + '. Remove those statevars and all their dependencies, or ' +
                         'supply a spline that includes ' + str(cutoff))
    # Make sure that spline has 3 dims if tdvs using concentration or f are requested
    reqF = [t for t in tdvSpec if t.reqF]
    reqM = [t for t in tdvSpec if t.reqM]
    if dimCt == 2 and (reqM or reqF):
        raise ValueError('You cannot calculate ' + pformat(set([t.name for t in (reqF + reqM)])) + ' with a spline ' +
                         'that does not include concentration. Remove those statevars and all their dependencies, ' +
                         'or supply a spline that includes concentration.')
    # Make sure that solvent molecular weight is provided if any tdvs that require MWv are requested
    reqMWv = [t for t in tdvSpec if t.reqMWv]
    if reqMWv and 'MW' not in gibbsSp:
        log.warning('When calculating ' + pformat([t.name for t in reqMWv]) +
                         ', solvent molecular weight defaults to the molecular weight of water.')
    # Make sure that solute molecular weight is provided if any tdvs that require MWu or f are requested
    reqMWu = [t for t in tdvSpec if t.reqMWu]
    MWu = None if 'MW' not in gibbsSp else np.array(gibbsSp['MW']) if isinstance(gibbsSp['MW'], float) else gibbsSp['MW']
    if (MWu is None or MWu.size != 2 or MWu[1] == 0) and (reqMWu or reqF): # currently only supporting single solute solutions
        raise ValueError('You cannot calculate ' + pformat([t.name for t in reqMWu] + [t.name for t in reqMWu]) + ' without ' +
                         'providing solute molecular weight.  Remove those statevars and all their dependencies, or ' +
                         'provide a valid value for the MWu parameter.')
    reqNu = [t for t in tdvSpec if t.reqNu]
    nu = None if 'nu' not in gibbsSp else gibbsSp['nu']
    if (nu == 0 or not nu) and reqNu:  # currently only supporting single solute solutions
        raise ValueError('You cannot calculate ' + pformat([t.name for t in reqNu]) + ' without ' +
                         'providing the number of ions in solution. ' + 'Remove those statevars and their dependencies, '
                         + 'or provide a valid value for the nu parameter.')
    # Make sure that all the PTM values fall inside the knot ranges of the spline
    # For single point, PTM is a tuple - just report its min and max as the single value
    # For a grid, PTM is an ndarray of sorted ndarrays, and so the min and max are the first and last elements
    # TODO: make this less explictly dependent on data types
    ptmranges = [(d[0], d[-1]) if isinstance(d, Iterable) else (d, d) for d in PTM]
    hasValsOutsideKnotRange = lambda kr, dr: dr[0] < kr[0] or dr[1] > kr[1]
    extrapolationDims = [i for i in range(0, dimCt) if hasValsOutsideKnotRange(knotranges[i], ptmranges[i])]
    if extrapolationDims:
        msg = ' '.join(
            ['Dimensions', pformat({'P' if d == iP else ('T' if d == iT else 'M') for d in extrapolationDims}),
             'contain values that fall outside the knot sequence for the given spline.'])
        if failOnExtrapolate:
            raise ValueError(msg)
        else:
            log.warning(msg)
    # Warn the user if the calculation results will take more than some factor times total virtual memory
    # TODO: make this less explicitly dependent on data types
    ptct = np.prod([len(d) if isinstance(d, Iterable) else 1 for d in PTM])  # the total number of points
    # For each point, the output will include 1 value for the PTM point itself plus 1 for each tdv
    outputSize = (len(tdvSpec) + 1) * ptct * floatSizeBytes
    if outputSize > virtual_memory().total * vmWarningFactor:
        log.warning(f'The projected output is more than {vmWarningFactor} times the total virtual memory for this machine.')
    return


def createThermodynamicStatesObj(tdvSpec, PTM, initializetdvs=False):
    """

    :param tdvSpec:         The tdvSpec indicating which tdvs will need to be stored in the object
    :param PTM:             The original PTM data, which will be stored in the object
    :param initializetdvs:  True if creating the ThermodynamicStates object for scatter data
                            Initializes output in the object by creating empty arrays matching the size of PTM
    :return:                A ThermodynamicStates object ready to store calculated thermodynamic state variables
    """
    flds = {t.name for t in tdvSpec} | {'PTM'}

    # Copy PTM so if you add a 0 concentration, you affect only the version in the output var
    # so later you can compare the original PTM to the one in tdvout.PTM to see if you need to remove the 0 M
    TDS = type('ThermodynamicStates', (object,),
               {fld: (np.copy(PTM) if fld == 'PTM'
                      else np.full((PTM.size,), np.nan, float) if initializetdvs
                      else None) for fld in flds})
    out = TDS()
    return out


def _ptmTuple2NestedArrays(PTM, dimCt, needs0M):
    """ Converts a PTM tuple to numpy ndarrays

    :return: PTM as nested arrays, possibly with introduced 0 M
    """
    out = np.empty(dimCt, object)
    for i in np.arange(0, dimCt):
        out[i] = np.empty((1,), float)
        # If a 0 concentration is used in PTM, replace with eps.
        out[i][0] = eps if i == iM and PTM[i] == 0 else PTM[i]
    return out


def _needsCutoff(dimCt, tdvSpec):
    """ Determines whether any tdvs require 0 concentration to be properly evaluated
    """
    return dimCt == 3 and any([t.reqCutoff for t in tdvSpec])

def _getPTMWithCutoff(PTM, cutoff):
    """ if all three dimensions happen to be the same size, numpy may create PTM as a 2D array
        rather than an array of three arrays, which is what we need.  The only place in the code
        this matters is when the size of one of the dimensions of PTM changes (as here).
    """
    newPTM = np.empty((3,), dtype=object)
    newPTM[iM] = np.insert(PTM[iM], 0, cutoff)
    newPTM[iP] = PTM[iP]
    newPTM[iT] = PTM[iT]
    return newPTM

def _removeCutoff(tdvout):
    """ If a "cutoff" concentration was added, take the values out. A cutoff concentration will always be prepended

    NOTE: This method changes the value of tdvout without returning it.

    :param tdvout:      The object with calculated thermodynamic properties
                        tdvout.PTM is the original PTM input provided by the caller
    """
    # Go through all calculated values and remove the first item from the M dimension
    # TODO: figure out why PTM doesn't show up in vars
    for p, v in vars(tdvout).items():
        if p != 'PTM':
            slc = [slice(None)] * len(v.shape)
            slc[iM] = slice(1, None)
            setattr(tdvout, p, v[tuple(slc)])


def _createGibbsDerivativesClass(tdvSpec):
    flds = {d for t in tdvSpec if t.reqDerivs for d in t.reqDerivs}
    return type('GibbsDerivatives', (object,), {d: None for t in tdvSpec if t.reqDerivs for d in t.reqDerivs})


def _buildDerivDirective(derivSpec, dimCt):
    """ Gets a list of the derivatives for relevant dimensions

    :param derivSpec: Derivatives specified
    :param dimCt: int. Dimension count (2 or 3)
    :return:
        out. Derivative specification
    """
    out = [statevars.defDer] * dimCt
    if derivSpec.wrtP: out[iP] = derivSpec.wrtP
    if derivSpec.wrtT: out[iT] = derivSpec.wrtT
    if derivSpec.wrtM: out[iM] = derivSpec.wrtM
    return out


def getDerivatives(gibbsSp, PTM, dimCt, tdvSpec, allowExtrapolations, verbose=False):
    """

    :param gibbsSp:     The Gibbs energy spline
    :param PTM:         The pressure, temperature[, molality] points at which the spine should be evaluated
    :param dimCt:       2 if a PT spline, 3 if a PTM spline
    :param tdvSpec:     The expanded TDVSpec describing the thermodynamic variables to be calculated
    :param allowExtrapolations False to return numpy.nan for all PTM values that fall outside the range of the spline
    :param verbose:     True to output timings
    :return:            Gibbs energy derivatives necessary to calculate the tdvs listed in the tdvSpec
    """
    GibbsDerivs = _createGibbsDerivativesClass(tdvSpec)
    out = GibbsDerivs()
    reqderivs = {d for t in tdvSpec for d in t.reqDerivs}  # Get set of derivative names that are needed
    getDerivSpec = lambda dn: next(d for d in statevars.derivatives if d.name == dn)
    for rd in reqderivs:
        derivDirective = _buildDerivDirective(getDerivSpec(rd), dimCt)
        start = time()
        setattr(out, rd, evalMultivarSpline(gibbsSp, PTM, derivDirective, allowExtrapolations))
        end = time()
        if verbose: _printTiming('deriv ' + rd, start, end)
    return out


def _getGriddedPTM(tdvSpec, PTM, verbose=False):
    if any([t.reqGrid for t in tdvSpec]):
        start = time()
        out = np.meshgrid(*PTM.tolist(), indexing='ij')  # Grid the dimensions of PTM
        end = time()
        if verbose: _printTiming('grid', start, end)
    else:
        out = []
    return out


def _getVolSolventInVolSlnConversion(MWu, PTM):
    """ Dimensionless conversion factor for how much of the volume of 1 kg of solution is really just the solvent
    :return: f
    """
    return 1 + MWu * PTM[iM]

def _isPointExtrapolation(dimCt, knots, PTM):
    for di in np.arange(dimCt):
        extrap = _isExtrapolation(knots[di], PTM[di])
        if bool(extrap.size):
            return True
    return False






def _printTiming(calcdesc, start, end):
    endDT = datetime.fromtimestamp(end)
    print(endDT.strftime('%H:%M:%S.%f'), ':\t', calcdesc, 'took', str(end - start), 'seconds to calculate')


#########################################
## Constants
#########################################
vmWarningFactor = 2  # Warn the user when size of output would exceed vmWarningFactor times total virtual memory
floatSizeBytes = int(np.finfo(float).bits / 8)
H20_MWv = 18.01528e-3