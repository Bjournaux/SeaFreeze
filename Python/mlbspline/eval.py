from scipy.interpolate import splev
import numpy as np

iT = 1

def evalMultivarSpline(spd, x, der=[]):
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
    :return:    a Numpy n-D array. If x is an ndarray, then output.shape == x.shape.
                If x is a tuple, then output.shape == tuple(np.ones(len(x), np.int)).
    """
    dimCt = spd['number'].size          # size of dim coefs
    if len(der) == 0:
        der = [0] * dimCt # default to 0th derivative for all dimensions

    if len(x) != dimCt:     # use len(x) because it works with both tuples and ndarrays
        raise ValueError("The dimensions of the evaluation points do not match the dimensions of the spline.")
    if len(der) != dimCt:
        raise ValueError("The dimensions of the derivative directives do not match the dimensions of the spline.")
    if not all((isinstance(i, int) and i >= 0) for i in der):
        raise ValueError('At least one derivative directive is not a non-negative integer.')
    # Use chain rule to get multiplicative factor for non-dimensional temperatures in B spline
    ndT1, ndT2 = (1, 1)
    if der[iT] > 0 and spd['ndT']:
        T_K = np.exp(x[iT])*spd['Tc']
        ndT1 = 1/T_K
        if der[iT] > 1:
            ndT2 = -1/T_K**2

    y = spd['coefs']
    # start at the last index and work downward to the first
    for di in range(dimCt - 1, -1, -1):
        xi = x[di]  # get x values for dimension being evaluated
        # wrap xi if necessary
        if not isinstance(xi, np.ndarray):
            xi = np.asarray([xi])
        tck = _getNextSpline(di, dimCt, spd, y)
        if di == iT and spd['ndT']:
            if der[iT] == 2 and spd['ndT']:
                # Use product rule to combine derivates in terms of non-dimensional T
                y = np.array(splev(xi, tck, der=2))*ndT1**2 + np.array(splev(xi, tck, der=1))*ndT2
            else:
                y = np.array(splev(xi, tck, der=der[di]))*ndT1
        else:
            y = np.array(splev(xi, tck, der=der[di]))
    # need to rearrange back to original order and shape
    return _getNextSpline(-1, dimCt, spd, y)[1]


def _getNextSpline(dimIdx, dimCt, spd, coefs):
    """ Get the tck (knots, coefs, and degree) for the next spline to be evaluated

    :param dimIdx:  The zero-based index of the dimension being evaluated
    :param dimCt:   The number of dimensions
    :param spd:     The spline being evaluated
    :param coefs:   The coefficients for this dimension - for this outermost dim, this is just the coefs from the spline
                    Otherwise, this is the output from the last call to interp.splev
    :return:        The [t,c,k] spline representation for the next dimension (see splev documentation)
    """
    li = dimCt - 1
    if li != dimIdx:
        coefs = np.moveaxis(coefs, li, 0)
    t = spd['knots'][dimIdx]
    k = spd['order'][dimIdx] - 1    # scipy wants the degree, not the order (MatLab gives the orders)
    return [t, coefs, k]
