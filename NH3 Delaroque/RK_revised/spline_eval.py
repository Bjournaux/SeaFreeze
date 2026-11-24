"""
Tensor-product B-spline evaluation with scipy.interpolate.BSpline.

SeaFreeze's mlbspline.eval uses scipy's `splev`, whose FIRST derivative is wrong for
splines of degree above 5 (order above 6): FITPACK only supports degree <= 5.  Values
and second derivatives are unaffected.  The earlier NH3-H2O spline "V1" is order 8 in
every dimension, so its density, sound speed and heat capacity must be evaluated with
this function, not with mlbspline.eval.evalMultivarSpline.
"""
import numpy as np
from scipy.interpolate import BSpline


def eval_bspline(sp, x, der):
    """Evaluate spline dict `sp` on the grid x = [x0, x1, ...] with derivative orders `der`."""
    c = np.asarray(sp["coefs"], float)
    for ax in range(len(x) - 1, -1, -1):
        b = BSpline(np.asarray(sp["knots"][ax], float), np.moveaxis(c, ax, 0), int(sp["order"][ax]) - 1,
                    extrapolate=False)
        c = np.moveaxis(b(np.asarray(x[ax], float), int(der[ax])), 0, ax)
    return c
