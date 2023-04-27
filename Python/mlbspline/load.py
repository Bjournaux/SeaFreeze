from scipy.io import loadmat
import numpy as np
import logging

log = logging.getLogger('mlbspline')
stream = logging.StreamHandler()
stream.setFormatter(logging.Formatter('[MLBspline %(levelname)s] %(message)s'))

def loadSpline(splineFile, splineVar=None):
    """Loads a spline from .mat format file
    Only MatLab files saved as either v7 or v7.3 are supported

    :param splineFile: full or relative path to Matlab file
    :param splineVar: variable to load from splineFile
    :return: a dict with the spline representation
    """
    raw = _getRaw(splineFile, splineVar)
    spd = getSplineDict(_stripNestingToFields(raw))
    # Knots can sometimes be imported as object type, which cannot be automatically cast to float
    if not np.issubdtype(spd['knots'][0].dtype, np.number) and not isinstance(spd['knots'][0][0], np.ndarray):
        spd['knots'] = spd['knots'].astype(np.float_)
    validateSpline(spd)
    return spd


def _getRaw(splineFile, splineVar=None):
    f = _loadFile(splineFile, splineVar)
    splineVar = _getCheckVar(f, splineVar)
    return f[splineVar]


def _loadFile(f, v=None):
    return loadmat(f, variable_names=None if v is None else [v], chars_as_strings=True)


def _getCheckVar(rawf, v=None):
    contents = [k for k in rawf.keys() if not k.startswith('__')]
    if v is None:
        if len(contents) == 1:
            v = contents[0]
        else:
            raise ValueError('The splineFile contains multiple variables: ' + ', '.join(contents) +
                             'Please provide the appropriate splineVar.')
    elif v not in contents:
        raise ValueError('The specified splineVar cannot be found in the specified file.  It only contains ' +
                         'variables ' + ', '.join(contents))
    return v


def _stripNestingToFields(src):
    lastFields = src.dtype.names
    while src[0].dtype.names == lastFields:
        src = src[0]
        lastFields = src.dtype.names
    return src


def _stripNestingToValue(src):
    # Numpy already automatically squeezes
    while isinstance(src, np.ndarray) and src.shape[0] == 1:
        src = src[0]
    return src


def getSplineDict(src):
    # All splines must contain the following fields
    out = {
        'form':     _stripNestingToValue(src['form']),
        'knots':    np.array([_stripNestingToValue(kd) for kd in _stripNestingToValue(src['knots'])], dtype=object),
        'number':   _stripNestingToValue(src['number']).astype(int),
        'order':    _stripNestingToValue(src['order']).astype(int),
        'dim':      _stripNestingToValue(src['dim']).astype(int),
        'coefs':    _stripNestingToValue(src['coefs']),
        'ndT':      False
    }

    # If dimensionless temperature is used, include the critical temperature in field ['Tc'] in the outputs to flag it
    if 'Tc' in src.dtype.names:
        out['Tc'] = _stripNestingToValue(src['Tc'])
        out['ndT'] = True

    # If number is a scalar, this is a 1D spline and some stuff needs to be re-wrapped for later code to work
    if not isinstance(out['number'], np.ndarray):
        out['number'] = np.array([out['number']])
        out['order'] = np.array([out['order']])
        # We can't just reshape or nest in one step, sadly. Numpy tries to be smart and creates a 2D array.
        # See https://stackoverflow.com/questions/7667799/numpy-object-arrays
        knots = np.empty(1, object)
        knots[0] = out['knots']
        out['knots'] = knots

    return out


def validateSpline(spd):
    """Checks the spline representation to make sure it has all the expected fields for evaluation
    Throws an error if anything is missing (one at a time)

    :param spd: a dict (output from getSplineDict) representing the spline
    """
    # Currently only supports b-splines
    if spd['form'] != 'B-':
        raise ValueError('These functions currently only support b-splines.')
    # Currently only supports 1-D output
    if spd['dim'] != 1:
        raise ValueError('These functions currently only support 1-D values.')
    # Need to have same size for knots, number, order
    if spd['number'].shape[0] != spd['knots'].shape[0] or spd['number'].shape != spd['order'].shape:
        raise ValueError('The spline''s knots, number, and order are inconsistent.')
    # Number of dims in coefs should be same as size of number (also makes it match knots and order b/c of prev check)
    if spd['number'].size != spd['coefs'].ndim:
        raise ValueError('The spline''s coefficients are inconsistent with its number.')
    # Each dim of coefs should have as many members as indicated by number
    if not (np.array(spd['coefs'].shape) == spd['number']).all():
        raise ValueError('At least one of the spline''s coefficients doesn''t match the corresponding number.')
    # Each dim of knots should have be of length dictated by corresponding number + order
    if not np.array(spd['number'] + spd['order'] == [k.size for k in spd['knots']]).all():
        raise ValueError('The number of knots does not correspond to the number and order for at least one variable.')
    return
