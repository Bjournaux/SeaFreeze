from datetime import datetime
from pprint import pformat
from time import time
from warnings import warn
from collections.abc import Iterable

import numpy as np
from psutil import virtual_memory

from mlbspline.eval import evalMultivarSpline
from lbftd import statevars
from lbftd.statevars import iT, iP, iM


def evalSolutionGibbsGrid(gibbsSp, PTM, *tdvSpec, MWv=18.01528e-3, MWu=None, failOnExtrapolate=True, verbose=False):
    """ Calculates thermodynamic variables for solutions based on a spline giving Gibbs energy
    This currently only supports single-solute solutions.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Warning: units must be as specified here because some conversions are hardcoded into this function.
        With the exception of pressure, units are SI.  Pressure is in MPa rather than Pa.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    :param gibbsSp: A B-spline  (in format given by loadSpline.loadMatSpline) for giving Gibbs energy (J/mg)
                    with dimensions pressure (MPa), temperature (K), and (optionally) molality (mol/kg), in the
                    order defined by statevars.iP, iT, and iM.
                    If molality is not provided, this function assumes that it is calculating
                    thermodynamic properties for a pure substance.
    :param PTM:     a Numpy ndarray of ndarrays with the conditions at which gibbsSp should be evaluated
                    the number of dimensions must be same as in the spline (PTM.size == gibbsSp['number'].size)
                    each of the inner ndarrays represents one of pressure (P), temperature (T), or molality (M)
                    and must be in the same order and units described in the notes for the gibbsSp parameter
                    Additionally, each dimension must be sorted from low to high values.
    :param MWv:     float with molecular weight of solvent (kg/mol).
                    Defaults to molecular weight of water (7 sig figs)
    :param MWu:     float with molecular weight of solute (kg/mol).
    :param failOnExtrapolate:   True if you want an error to appear if PTM includes values that fall outside the knot
                    sequence of gibbsSp.  If False, throws a warning rather than an error, and
                    proceeds with the calculation.
    :param verbose: boolean indicating whether to print status updates, warnings, etc.
    :param tdvSpec: iterable indicating the thermodynamic variables to be calculated
                    elements can be either strings showing the names (full list at statevars.statevarnames)
                    or TDV objects from statevars.statevars.
                    If not provided, this function will calculate the variables in statevars.tdvsPTOnly for a PT spline,
                    and those in statevars.statevars for a PTM spline.
                    Any other args provided will result in an error.
                    See the README for units.
    :return:        a named tuple with the requested thermodynamic variables as named properties
                    matching the statevars requested in the *tdvSpec parameter of this function
                    the output will also include P, T, and M (if provided) properties
    """
    [dimCt, tdvSpec] = _parseInput(gibbsSp, *tdvSpec)
    _checkInputs(gibbsSp, dimCt, tdvSpec, PTM, MWv, MWu, failOnExtrapolate)
    tdvout = createThermodynamicStatesObj(tdvSpec, PTM)  # this has the original PTM
    if _needs0M(dimCt, tdvSpec) and PTM[iM][0] != 0:
        PTM[iM] = np.insert(PTM[iM], 0, 0)
    tdvout = _evalInternal(gibbsSp, tdvSpec, PTM, MWu, MWv, tdvout, dimCt, verbose)
    return tdvout


def evalSolutionGibbsScatter(gibbsSp, PTM, *tdvSpec, MWv=18.01528e-3, MWu=None, failOnExtrapolate=True, verbose=False):
    """ Calculates thermodynamic variables for solutions based on a spline giving Gibbs energy
    This currently only supports single-solute solutions.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Warning: units must be as specified here because some conversions are hardcoded into this function.
        With the exception of pressure, units are SI.  Pressure is in MPa rather than Pa.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    :param gibbsSp: A B-spline  (in format given by loadSpline.loadMatSpline) for giving Gibbs energy (J/mg)
                    with dimensions pressure (MPa), temperature (K), and (optionally) molality (mol/kg), in the
                    order defined by statevars.iP, iT, and iM.
                    If molality is not provided, this function assumes that it is calculating
                    thermodynamic properties for a pure substance.
    :param PTM:     a Numpy ndarray of tuples with the points at which gibbsSp should be evaluated
                    the number of elements in each tuple must be same as in the spline (PTM[i].size == gibbsSp['number'].size)
                    each tuples gives the pressure (P), temperature (T), and (optionally) molality (M)
                    and must be in the same order and units described in the notes for the gibbsSp parameter
                    No sorting is required.
    :param MWv:     float with molecular weight of solvent (kg/mol).
                    Defaults to molecular weight of water (7 sig figs)
    :param MWu:     float with molecular weight of solute (kg/mol).
    :param failOnExtrapolate:   True if you want an error to appear if PTM includes values that fall outside the knot
                    sequence of gibbsSp.  If False, throws a warning rather than an error, and
                    proceeds with the calculation.
    :param verbose: boolean indicating whether to print status updates, warnings, etc.
    :param tdvSpec: iterable indicating the thermodynamic variables to be calculated
                    elements can be either strings showing the names (full list at statevars.statevarnames)
                    or TDV objects from statevars.statevars.
                    If not provided, this function will calculate the variables in statevars.tdvsPTOnly for a PT spline,
                    and those in statevars.statevars for a PTM spline.
                    Any other args provided will result in an error.
                    See the README for units.
    :return:        a named tuple with the requested thermodynamic variables as named properties
                    matching the statevars requested in the *tdvSpec parameter of this function
                    Each statevar property will contain an ndarray of values corresponding to the points listed in PTM.
                    The output also contains a PTM property for convenience.
    """
    [dimCt, tdvSpec] = _parseInput(gibbsSp, *tdvSpec)
    needs0M = _needs0M(dimCt, tdvSpec)
    fakePTM = _makeFakePTMGrid(dimCt, PTM)  # fakes grid-style input to do all input checks before processing
    _checkInputs(gibbsSp, dimCt, tdvSpec, fakePTM, MWv, MWu, failOnExtrapolate)

    tdvout = createThermodynamicStatesObj(tdvSpec, PTM, initializetdvs=True)
    # point by point
    for p in np.arange(0, PTM.size):
        wPTM = _ptmTuple2NestedArrays(PTM[p], dimCt, needs0M)
        tempout = createThermodynamicStatesObj(tdvSpec, wPTM)
        # prepend a 0 concentration if one is needed by any of the quantities being calculated
        if needs0M and wPTM[iM] != 0:
            wPTM[iM] = np.insert(wPTM[iM], 0, 0)
        _evalInternal(gibbsSp, tdvSpec, wPTM, MWu, MWv, tempout, dimCt, verbose=False)
        # copy temp output to the real output
        for t in tdvSpec:
            tolastpt = getattr(tdvout, t.name)
            tolastpt[p] = getattr(tempout, t.name)
            setattr(tdvout, t.name, tolastpt)
    return tdvout


def _evalInternal(gibbsSp, tdvSpec, PTM, MWu, MWv, tdvout, dimCt, verbose=False):
    derivs = getDerivatives(gibbsSp, PTM, dimCt, tdvSpec, verbose)
    gPTM = _getGriddedPTM(tdvSpec, PTM, verbose) if any([tdv.reqGrid for tdv in tdvSpec]) else None
    f = _getVolSolventInVolSlnConversion(MWu, PTM) if any([tdv.reqF for tdv in tdvSpec]) else None

    # calculate thermodynamic variables and store in appropriate fields in tdvout
    completedTDVs = set()  # list of completed tdvs
    while len(completedTDVs) < len(tdvSpec):
        # get tdvs that either have empty reqTDV or all of those tdvs have already been calculated
        # but that are not in completedTDVs themselves (don't recalculate anything)
        tdvsToEval = tuple(t for t in tdvSpec if
                           t.name not in completedTDVs and (not t.reqTDV or not t.reqTDV.difference(completedTDVs)))
        for t in tdvsToEval:
            args = _buildEvalArgs(t, derivs, PTM, gPTM, MWv, MWu, tdvout, gibbsSp, f)
            start = time()
            setattr(tdvout, t.name, t.calcFn(**args))  # calculate the value and set it in the output
            end = time()
            if verbose: _printTiming('tdv ' + t.name, start, end)
            completedTDVs.add(t.name)
    _remove0M(tdvout, PTM)
    return tdvout


def _parseInput(gibbsSp, *tdvSpec):
    dimCt = gibbsSp['number'].size
    # expand spec to add dependencies (or set to default spec if no spec given)
    origSpec = tdvSpec
    tdvSpec = statevars.expandTDVSpec(tdvSpec, dimCt)
    addedTDVs = [s.name for s in tdvSpec if s.name not in origSpec]
    if origSpec and addedTDVs:  # the original spec was not empty and more tdvs were added
        print('NOTE: The requested thermodynamic variables depend on the following variables, which will be ' +
              'included as properties of the output object: ' + pformat(addedTDVs))
    return dimCt, tdvSpec


def _buildEvalArgs(tdv, derivs, PTM, gPTM, MWv, MWu, tdvout, gibbsSp, f):
    args = dict()
    if tdv.reqDerivs: args[tdv.parmderivs] = derivs
    if tdv.reqGrid:   args[tdv.parmgrid] = gPTM
    if tdv.reqMWv:    args[tdv.parmMWv] = MWv
    if tdv.reqMWu:    args[tdv.parmMWu] = MWu
    if tdv.reqTDV:    args[tdv.parmtdv] = tdvout
    if tdv.reqSpline: args[tdv.parmspline] = gibbsSp
    if tdv.reqPTM:    args[tdv.parmptm] = PTM
    if tdv.reqF:      args[tdv.parmf] = f
    return args


def _makeFakePTMGrid(dimCt, PTM):
    out = np.empty(dimCt, np.object)
    for d in np.arange(0, dimCt):
        u = {t[d] for t in PTM}                 # get unique values in each dimension
        out[d] = np.array(sorted(u))            # sort the values in ascending order
    return out


def _checkInputs(gibbsSp, dimCt, tdvSpec, PTM, MWv, MWu, failOnExtrapolate):
    """ Checks error conditions before performing any calculations, throwing an error if anything doesn't match up
      - Ensures necessary data is available for requested statevars (check req parameters for statevars)
      - Throws a warning (or error if failOnExtrapolate=True) if spline is being evaluated on values
        outside the range of its knots
      - Warns the user if the requested output would exceed vmWarningFactor times the amount of virtual memory
      Note that the mlbspline.eval module performs additional checks (like dimension count mismatches
    """
    knotranges = [(k[0], k[-1]) for k in gibbsSp['knots']]
    # if a PTM spline, make sure the spline's concentration dimension starts with 0 if any tdv has req0M=True
    # Note that the evalGibbs functions elsewhere handle the case where PTM does not include 0 concentration.
    req0M = [t for t in tdvSpec if t.req0M]
    if dimCt == 3 and req0M and knotranges[iM][0] != 0:
            raise ValueError('You cannot calculate ' + pformat([t.name for t in req0M]) + ' with a spline that does ' +
                             'not include 0 concentration. Remove those statevars and all their dependencies, or ' +
                             'supply a spline that includes 0 concentration.')
    # make sure that spline has 3 dims if tdvs using concentration or f are requested
    reqF = [t for t in tdvSpec if t.reqF]
    reqM = [t for t in tdvSpec if t.reqM]
    if dimCt == 2 and (reqM or reqF):
        raise ValueError('You cannot calculate ' + pformat(set([t.name for t in (reqF + reqM)])) + ' with a spline ' +
                         'that does not include concentration. Remove those statevars and all their dependencies, ' +
                         'or supply a spline that includes concentration.')
    # make sure that solvent molecular weight is provided if any tdvs that require MWv are requested
    reqMWv = [t for t in tdvSpec if t.reqMWv]
    if (MWv == 0 or not MWv) and reqMWv:
        raise ValueError('You cannot calculate ' + pformat([t.name for t in reqMWv]) + ' without ' +
                         'providing solvent molecular weight.  Remove those statevars and all their dependencies, or ' +
                         'provide a valid value for the MWv parameter.')
    # make sure that solute molecular weight is provided if any tdvs that require MWu or f are requested
    reqMWu = [t for t in tdvSpec if t.reqMWu]
    if (MWu == 0 or not MWu) and (reqMWu or reqF):
        raise ValueError('You cannot calculate ' + pformat([t.name for t in set(reqMWu + reqF)]) + ' without ' +
                         'providing solute molecular weight.  Remove those statevars and all their dependencies, or ' +
                         'provide a valid value for the MWu parameter.')
    # make sure that all the PTM values fall inside the knot ranges of the spline
    # for single point, PTM is a tuple - just report its min and max as the single value
    # for a grid, PTM is an ndarray of sorted ndarrays, and so the min and max are the first and last elements
    # TODO: make this less explictly dependent on data types
    ptmranges = [(d[0], d[-1]) if isinstance(d, Iterable) else (d, d) for d in PTM]
    hasValsOutsideKnotRange = lambda kr, dr: dr[0] < kr[0] or dr[1] > kr[1]
    extrapolationDims = [i for i in range(0, dimCt) if hasValsOutsideKnotRange(knotranges[i], ptmranges[i])]
    if extrapolationDims:
        msg = ' '.join(['Dimensions',pformat({'P' if d == iP else ('T' if d == iT else 'M') for d in extrapolationDims}),
                'contain values that fall outside the knot sequence for the given spline,',
                'which will result in extrapolation, which may not produce meaningful values.'])
        if failOnExtrapolate:
            raise ValueError(msg)
        else:
            warn(msg)
    # warn the user if the calculation results will take more than some factor times total virtual memory
    # TODO: make this less explicitly dependent on data types
    ptct = np.prod([len(d) if isinstance(d, Iterable) else 1 for d in PTM])       # the total number of points
    # for each point, the output will include 1 value for the PTM point itself plus 1 for each tdv
    outputSize = (len(tdvSpec) + 1) * ptct * floatSizeBytes
    if outputSize > virtual_memory().total * vmWarningFactor:
        warn('The projected output is more than {0} times the total virtual memory for this machine.'.format(vmWarningFactor))
    return


def createThermodynamicStatesObj(tdvSpec, PTM, initializetdvs=False):
    """

    :param tdvSpec:         the tdvSpec indicating which tdvs will need to be stored in the object
    :param PTM:             the original PTM data, which will be stored in the object
    :param initializetdvs:  True if creating the ThermodynamicStates object for scatter data
                            initializes output in the object by creating empty arrays matching the size of PTM
    :return:                a ThermodynamicStates object ready to store calculated thermodynamic state variables
    """
    flds = {t.name for t in tdvSpec} | {'PTM'}

    # copy PTM so if you add a 0 concentration, you affect only the version in the output var
    # so later you can compare the original PTM to the one in tdvout.PTM to see if you need to remove the 0 M
    TDS = type('ThermodynamicStates', (object,),
               {fld: (np.copy(PTM) if fld == 'PTM' else np.empty((PTM.size,), float) if initializetdvs else None) for fld in flds})
    out = TDS()
    return out


def _ptmTuple2NestedArrays(PTM, dimCt, needs0M):
    """ Converts a PTM tuple to numpy ndarrays

    :return:  PTM as nested arrays, possibly with introduced 0 M
    """
    out = np.empty(dimCt, np.object)
    for i in np.arange(0, dimCt):
        out[i] = np.empty((1,), float)
        out[i][0] = PTM[i]
    return out


def _needs0M(dimCt, tdvSpec):
    """ Determines whether any tdvs require 0 concentration to be properly evaluated
    """
    return dimCt == 3 and any([t.req0M for t in tdvSpec])


def _remove0M(tdvout, wrappedPTM):
    """ If a 0 concentration was added, take it back out.
    NOTE: This method changes the value of tdvout without returning it.

    :param tdvout:      The object with calculated thermodynamic properties
                        tdvout.PTM is the original PTM input provided by the caller
    :param wrappedPTM:  The wrapped PTM which may include
    """
    if wrappedPTM.size == 3 and wrappedPTM[iM][0] == 0:
        try:    # tdvout.PTM is either a tuple or another set of nested arrays
            firstM = tdvout.PTM[iM][0]
        except IndexError:
            firstM = tdvout.PTM[iM]
        if firstM != 0:  # this means the wrapped PTM had a 0 concentration prepended
            # go through all calculated values and remove the first item from the M dimension
            # TODO: figure out why PTM doesn't show up in vars
            for p, v in vars(tdvout).items():
                slc = [slice(None)] * len(v.shape)
                slc[iM] = slice(1, None)
                setattr(tdvout, p, v[tuple(slc)])


def _createGibbsDerivativesClass(tdvSpec):
    flds = {d for t in tdvSpec if t.reqDerivs for d in t.reqDerivs}
    return type('GibbsDerivatives', (object,), {d: None for t in tdvSpec if t.reqDerivs for d in t.reqDerivs})


def _buildDerivDirective(derivSpec, dimCt):
    """ Gets a list of the derivatives for relevant dimensions
    """
    out = [statevars.defDer] * dimCt
    if derivSpec.wrtP: out[iP] = derivSpec.wrtP
    if derivSpec.wrtT: out[iT] = derivSpec.wrtT
    if derivSpec.wrtM: out[iM] = derivSpec.wrtM
    return out


def getDerivatives(gibbsSp, PTM, dimCt, tdvSpec, verbose=False):
    """

    :param gibbsSp:     The Gibbs energy spline
    :param PTM:         The pressure, temperature[, molality] points at which the spine should be evaluated
    :param dimCt:       2 if a PT spline, 3 if a PTM spline
    :param tdvSpec:     The expanded TDVSpec describing the thermodynamic variables to be calculated
    :param verbose:     True to output timings
    :return:            Gibbs energy derivatives necessary to calculate the tdvs listed in the tdvSpec
    """
    GibbsDerivs = _createGibbsDerivativesClass(tdvSpec)
    out = GibbsDerivs()
    reqderivs = {d for t in tdvSpec for d in t.reqDerivs}   # get set of derivative names that are needed
    getDerivSpec = lambda dn: next(d for d in statevars.derivatives if d.name == dn)
    for rd in reqderivs:
        derivDirective = _buildDerivDirective(getDerivSpec(rd), dimCt)
        start = time()
        setattr(out, rd, evalMultivarSpline(gibbsSp, PTM, derivDirective))
        end = time()
        if verbose: _printTiming('deriv '+rd, start, end)
    return out


def _getGriddedPTM(tdvSpec, PTM, verbose=False):
    if any([t.reqGrid for t in tdvSpec]):
        start = time()
        out = np.meshgrid(*PTM.tolist(), indexing='ij')    # grid the dimensions of PTM
        end = time()
        if verbose: _printTiming('grid', start, end)
    else:
        out = []
    return out


def _getVolSolventInVolSlnConversion(MWu, PTM):
    """ Dimensionless conversion factor for how much of the volume of 1 kg of solution is really just the solvent
    :return:    f
    """
    return 1 + MWu * PTM[iM]


def _printTiming(calcdesc, start, end):
    endDT = datetime.fromtimestamp(end)
    print(endDT.strftime('%H:%M:%S.%f'), ':\t', calcdesc,'took',str(end-start),'seconds to calculate')


#########################################
## Constants
#########################################
vmWarningFactor = 2         # warn the user when size of output would exceed vmWarningFactor times total virtual memory
floatSizeBytes = int(np.finfo(float).bits / 8)
















