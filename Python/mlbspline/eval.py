from scipy.interpolate import splev
import numpy as np
import logging

log = logging.getLogger('mlbspline')
stream = logging.StreamHandler()
stream.setFormatter(logging.Formatter('[MLBspline %(levelname)s] %(message)s'))

iT = 1

def evalMultivarSpline(spd, x, der=None, allowExtrapolations=True):
    """ Performs recursive evaluation of b-spline for the given independent values
    For now, assumes 1-D spline (y for each n-D x is scalar)
    x and spd['coefs'] must have the same number of dimensions

    :param spd: a B-spline definition as given by load.loadSpline
    :param x:   Either
                    - a Numpy n-D array with the points at which the spline should be evaluated
                    or
                    - a tuple representing a single point at which the spline should be evaluated
                Either way, the number of dimensions must be same as in the spline (i.e.,
                len(x) == spd['number'].size).  Dimensions must be in the same order as in the spline.
    :param der: (optional) a list of the derivative levels for evaluation, with each value representing
                the derivative for the dimension at the corresponding index in x
                if provided, values must be non-negative integers, and
                the number of dimensions must be the same as in the spline (len(der) == spd['number'].size)
    :param allowExtrapolations: if False, any values of x that fall outside the knot range of the spline
                given by spd will return numpy.nan.  Note that you can initially get extrapolations just to see what
                they look like and clean them up later using _setExtrapolationsToNan.
                if True (default), extrapolations will be calculated according to the
                spline's coefficients at its boundaries.  Not recommended, but used as default to preserve behavior.
    :return:    a Numpy n-D array. If x is a ndarray, then output.shape == x.shape.
                If x is a tuple, then output.shape == tuple(np.ones(len(x), np.int)).
    """
    if der is None:
        der = []

    dimCt = _getSplineDimCount(spd)
    if len(der) == 0:
        der = [0] * dimCt  # Default to 0th derivative for all dimensions

    if len(x) != dimCt:  # Use len(x) because it works with both tuples and ndarrays
        raise ValueError('The dimensions of the evaluation points do not match the dimensions of the spline.')
    if len(der) != dimCt:
        raise ValueError('The dimensions of the derivative directives do not match the dimensions of the spline.')
    if not all((isinstance(i, int) and i >= 0) for i in der):
        raise ValueError('At least one derivative directive is not a non-negative integer.')
    y = spd['coefs']
    # Start at the last index and work downward to the first
    for di in range(dimCt - 1, -1, -1):
        xi = x[di]  # Get x values for dimension being evaluated
        # Wrap xi if necessary
        if not isinstance(xi, np.ndarray):
            xi = np.asarray([xi])
        tck = _getNextSpline(di, dimCt, spd, y)
        # ext = 0 means extrapolations are calculated, where ext=1 means that values of x
        # that fall outside the knot range will return 0 as the value.  This is cleaned up later,
        # as we want to return numpy.nan instead, but some cycles may be saved, as presumably the
        # underlying code will not bother to perform the calculations (may want to test this
        # presumption at some point if this function requires optimization).
        y = np.array(splev(xi, tck, der=der[di], ext=0 if allowExtrapolations else 1))
    # Need to rearrange back to original order and shape
    out = _getNextSpline(-1, dimCt, spd, y)[1]
    if not allowExtrapolations:
        out = _setExtrapolationsToNan(spd, x, out)
    return out


def _getNextSpline(dimIdx, dimCt, spd, coefs):
    """ Get the tck (knots, coefs, and degree) for the next spline to be evaluated
    Note: For a multivariate spline, the outermost (highest index) dimension is evaluated first,
    and the output from that spline gives the coefficients for the next innermost spline.  Actual
    values for the overall spline are not given until the 0th dimension is evaluated.  This function
    gets the knots, coefficients, and degree for the next innermost spline.

    :param dimIdx:  The zero-based index of the dimension being evaluated
    :param dimCt:   The number of dimensions
    :param spd:     The spline being evaluated
    :param coefs:   The coefficients for this dimension - for this outermost dim, this is just the coefs from the spline
                    Otherwise, this is the output from the last call to interp.splev
    :return:        The [t,c,k] spline representation for the next dimension (see splev documentation),
                    where t represents the knots, c the coefficients, and k the order
    """
    li = dimCt - 1
    if li != dimIdx:
        coefs = np.moveaxis(coefs, li, 0)
    t = spd['knots'][dimIdx]
    k = spd['order'][dimIdx] - 1  # Scipy wants the degree, not the order (MatLab gives the orders)
    return [t, coefs, k]

def _isExtrapolation(knots, x):
    """ Gets indices for values that fall outside the knot range of one dimension of a
    multivariate spline

    :param knots:   The knots for the dimension of the spline being evaluated
    :param x:       Input values for the dimension as described in evalMultivarSpline
    """
    return np.argwhere(np.logical_or(x < knots.min(), x > knots.max())).squeeze()

def _setExtrapolationsToNan(spd, x, y):
    """ Returns a copy of y that has nan instead of the original value for all points where the x
     value falls outside the range of the spline.
    This function is independent of evalMultivarSpline() so could be used to cleanse data
    after processing if extrapolations were initially allowed.

    :param spd:      The spline that was evaluated with x input values to produce the output y
    :param x:       Input values as in evalMultivarSpline
    :param y:       Output values calculated by scipy for the provided values of x
    :return:        a copy of y with values changed to numpy.nan if the corresponding value of x falls outside the
                    knot range for that dimension of the spline.
    """
    dc = _getSplineDimCount(spd)
    out = y.copy()
    for di in range(0, dc): # do each dimension one at a time
        xi = x[di]
        ki = spd['knots'][di]
        si = [slice(None)] * dc
        si[di] = _isExtrapolation(ki, xi)
        out[tuple(si)] = np.nan
    return out


def _getSplineDimCount(sp):
    return sp['number'].size  # Size of dim coefs



