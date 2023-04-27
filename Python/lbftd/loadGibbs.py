from math import isnan
from numbers import Number
import logging

import numpy as np

from mlbspline import load

log = logging.getLogger('lbftd')
stream = logging.StreamHandler()
stream.setFormatter(logging.Formatter('[LBFTD %(levelname)s] %(message)s'))

def loadGibbsSpline(splineFile, splineVar=None):
    """Loads a Gibbs energy spline from .mat format file
    A Gibbs energy spline should be a Matlab struct that includes the following fields:
      - sp is the main spline itself.  It must be a 2D (PT) or 3D (PTX) spline as outlined in
            mlbspline.eval.evalMultivarSpline
      - MW is a list of the molecular weights of each species in the solution, measured in kg/mol,
            such that MW[0] is the molecular weight of the solvent (even if not used, i.e., if reqMWv is false)
            and MW[1] is the molecular weight of the solute
            Field must be present and can have length of 1 (pure substance) or 2 (single solute solution).
      - nu is the number of ions in solution for the solute.  All values must be positive integers.
            If MW.size == 1 (pure substance), the value will be ignored so the field can be absent --
            otherwise the field must be present and have a numeric value.

    :param splineFile:  full or relative path to Matlab file
    :param splineVar:   variable to load from splineFile.
                        If not provided, the splineFile must contain exactly one variable
    :return:            a dict with the Gibbs energy spline representation required by evalGibbs functions
    """
    raw = load._stripNestingToFields(load._getRaw(splineFile, splineVar))
    sp = load.getSplineDict(load._stripNestingToValue(raw['sp']))
    load.validateSpline(sp)
    MW = _getMW(raw)
    out = {
        'sp':   sp,
        'MW':   MW
    }
    if MW.size > 1:     #  only add nu if dealing with a solution rather than a pure substance
        out['nu'] = _getnu(raw)
    return out

def _getMW(raw):
    try:
        MW = load._stripNestingToValue(raw['MW'])
    except KeyError:
        log.warning('Could not load MW - defaulting to empty value - will not be able to calculate some thermodynamic values')
        MW = np.empty(0)
    MW = np.array([MW]) if isinstance(MW, Number) else MW   # wrap scalar value
    if MW.size > 2:
        raise ValueError('MW has too many elements, as multi-solute solutions are not currently supported.')
    if any([not isinstance(mw, Number) or isnan(mw) for mw in MW]):
        raise ValueError('At least one value of MW is not numeric, so the entire value is ignored.')
    return MW


def _getnu(raw):
    # currently only supporting a single solute solution so this is just an integer
    nu = load._stripNestingToValue(raw['nu'])
    if not isinstance(nu, Number) or nu != int(nu):
        raise ValueError('At least one value in nu is not an integer.')
    return int(nu)

