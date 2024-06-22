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

    :param splineFile:  Full or relative path to Matlab file
    :param splineVar:   Variable to load from splineFile.
                        If not provided, the splineFile must contain exactly one variable
    :return:            A dict with the Gibbs energy spline representation required by evalGibbs functions
    """
    raw = load._stripNestingToFields(load._getRaw(splineFile, splineVar))
    sp = load.getSplineDict(load._stripNestingToValue(raw['sp']))
    sp['MW'] = _getMW(sp)
    return {'sp': sp} # TODO: this is for legacy reasons, clean it up

def _getMW(src):
    try:
        MW = src['MW']
    except KeyError:
        log.warning('Could not load MW - defaulting to empty value - will not be able to calculate some thermodynamic values')
        MW = np.empty(0)
    MW = np.array([MW]) if isinstance(MW, Number) else MW  # Wrap scalar value
    if MW.size > 2:
        raise ValueError('MW has too many elements, as multi-solute solutions are not currently supported.')
    if any([not isinstance(mw, Number) or isnan(mw) for mw in MW]):
        raise ValueError('At least one value of MW is not numeric. The entire value will be ignored.')
    return MW


def _getnu(src):
    # Currently only supports a single solute solution so this is just an integer
    nu = src['nu']
    if not isinstance(nu, Number) or nu != int(nu):
        raise ValueError('At least one value in nu is not an integer.')
    return int(nu)

def _getcutoff(src):
    cutoff = src['cutoff']
    return float(cutoff)
