"""
Build a SeaFreeze Gibbs-energy spline from any function G(P, T[, m]).

This is the entry point for adding a new material (in particular a new aqueous
binary) to SeaFreeze:

    1. write G(P, T, m), the specific Gibbs energy of the solution
    2. sp = build_gibbs_spline(G, P_sites, T_sites, m_sites)
    3. save_gibbs_spline(sp, '<splines>/<folder>/<file>.mat', MW=..., nu=...)
    4. seafreeze.register_solution('<name>', '<solute>', MW, nu, cutoff,
                                   spline=('<folder>', '<file>.mat'))

Conventions (those of lbftd, which evaluates the spline):

    G   J per kg of *solution*
    P   MPa
    T   K
    m   mol of solute per kg of water (molality); omit for a pure substance

The result is a tensor-product B-spline in the same B-form as the splines shipped
with SeaFreeze, interpolating G exactly at the given sites.  All properties
(density, heat capacity, chemical potentials, ...) are derivatives of it, so the
sites must be dense enough for the *derivatives* that matter to be accurate; use
``check_gibbs_spline`` to measure that against the generating function.
"""

import os

import numpy as np
from scipy.interpolate import make_interp_spline
from scipy.io import savemat

DEFAULT_ORDER = (6, 6, 4)        # P, T, m : same as the SeaFreeze NaCl(aq) splines


# ---------------------------------------------------------------------------
# Sites and knots
# ---------------------------------------------------------------------------
def dilute_sites(m_max, n, power=2.0, m_min=0.0):
    """Molality sites clustered near the dilute end: m_min + (m_max-m_min) (i/(n-1))**power.

    The ideal-mixing term of a solution behaves like m ln m, whose slope diverges
    at m = 0; sites (hence knots) must be much denser there than at high molality.
    """
    u = np.linspace(0.0, 1.0, int(n))
    return m_min + (m_max - m_min) * u**power


def greville_sites(knots, order):
    """Greville abscissae of a knot vector: the sites at which interpolation in
    that spline space is always well posed (and reproduces a spline of the space)."""
    t = np.asarray(knots, float)
    k = int(order) - 1
    n = t.size - int(order)
    sites = np.array([t[i + 1:i + k + 1].mean() for i in range(n)])
    # the mean of repeated end knots can leave the interval by one rounding error
    return np.clip(sites, t[k], t[n])


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------
def build_gibbs_spline(G, P, T, m=None, order=None, knots=None):
    """Interpolating tensor-product B-spline of G on the given sites.

    :param G:      callable G(P, T) or G(P, T, m) -> J/kg.  It is called once with
                   broadcastable arrays of shape (nP,1[,1]), (1,nT[,1]), (1,1,nm) and
                   must return an array of shape (nP, nT[, nm]).  A precomputed
                   array of that shape is also accepted.
    :param P, T, m: 1-D strictly increasing data sites [MPa, K, mol/kg]; m=None for
                   a pure substance.
    :param order:  spline order (degree + 1) per dimension; default (6, 6, 4).
    :param knots:  optional list of full knot vectors, one per dimension (None entries
                   use the default not-a-knot vector of the sites).  With explicit
                   knots the number of sites must equal len(knots) - order.
    :return:       spline dict with the fields SeaFreeze expects:
                   form, knots, number, order, dim, coefs.
    """
    axes = [np.asarray(P, float).ravel(), np.asarray(T, float).ravel()]
    if m is not None:
        axes.append(np.asarray(m, float).ravel())
    nd = len(axes)
    order = tuple(int(o) for o in (order if order is not None else DEFAULT_ORDER[:nd]))
    if len(order) != nd:
        raise ValueError(f"order must have {nd} entries, got {len(order)}.")
    knots = list(knots) if knots is not None else [None] * nd
    if len(knots) != nd:
        raise ValueError(f"knots must have {nd} entries, got {len(knots)}.")
    for name, x, o in zip('PTm', axes, order):
        if x.size < o:
            raise ValueError(f"{name}: need at least {o} sites for order {o}, got {x.size}.")
        if np.any(np.diff(x) <= 0):
            raise ValueError(f"{name}: sites must be strictly increasing.")

    if callable(G):
        mesh = [x.reshape([-1 if i == j else 1 for j in range(nd)]) for i, x in enumerate(axes)]
        values = np.asarray(G(*mesh), float)
        values = np.broadcast_to(values, tuple(x.size for x in axes)).copy()
    else:
        values = np.array(G, float)
    if values.shape != tuple(x.size for x in axes):
        raise ValueError(f"G returned shape {values.shape}, expected {tuple(x.size for x in axes)}.")
    if not np.all(np.isfinite(values)):
        raise ValueError("G is not finite at every site; restrict the sites to the domain of G.")

    # Tensor-product interpolation: interpolate along one dimension at a time; the
    # coefficients of one pass are the data of the next.
    coefs = values
    out_knots = np.empty(nd, dtype=object)
    for ax, (x, o, t) in enumerate(zip(axes, order, knots)):
        kw = {} if t is None else {'t': np.asarray(t, float)}
        spl = make_interp_spline(x, coefs, k=o - 1, axis=ax, **kw)
        # scipy stores the coefficients with the interpolation axis first
        coefs = np.moveaxis(np.asarray(spl.c), 0, ax)
        out_knots[ax] = np.asarray(spl.t, float)

    return {
        'form': 'B-',
        'knots': out_knots,
        'number': np.array(coefs.shape, dtype=int),
        'order': np.array(order, dtype=int),
        'dim': 1,
        'coefs': coefs,
    }


# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
def save_gibbs_spline(sp, filename, MW=None, nu=None, cutoff=None, **metadata):
    """Write a spline to a .mat file readable by SeaFreeze (variable ``sp``).

    :param MW:       solute molar mass [kg/mol] (solutions); omit for a pure substance.
    :param nu:       ions per formula unit of solute (solutions).
    :param cutoff:   dilute molality cutoff used for apparent properties [mol/kg].
    :param metadata: extra fields stored with the spline (strings or numbers), e.g.
                     model='...', calibrated='muw, aw', source='...'.
    """
    nd = len(sp['knots'])
    cell = np.empty((1, nd), dtype=object)
    for i, k in enumerate(sp['knots']):
        cell[0, i] = np.asarray(k, float).reshape(1, -1)
    out = {
        'form': 'B-',
        'knots': cell,
        'number': np.asarray(sp['number'], float).reshape(1, -1),
        'order': np.asarray(sp['order'], float).reshape(1, -1),
        'dim': 1.0,
        'coefs': np.asarray(sp['coefs'], float),
    }
    if MW is not None:
        out['MW'] = float(MW)
    if nu is not None:
        out['nu'] = float(nu)
    if cutoff is not None:
        out['cutoff'] = float(cutoff)
    for key, val in metadata.items():
        if key in out:
            raise ValueError(f"metadata key {key!r} clashes with a spline field.")
        out[key] = val
    folder = os.path.dirname(os.path.abspath(filename))
    os.makedirs(folder, exist_ok=True)
    savemat(filename, {'sp': out}, do_compression=True)
    return filename


# ---------------------------------------------------------------------------
# Check against the generating function
# ---------------------------------------------------------------------------
def check_gibbs_spline(sp, G, P, T, m=None, MW=None):
    """Compare a spline with its generating function at test points (ideally *not*
    the building sites).  Derivatives of G are taken by central differences.

    :return: dict of maximum absolute errors
             'G' [J/kg], 'rho' [kg/m3] and 'Cp' [J/kg/K], and for a solution with MW
             given also 'muw' [J/mol of H2O].
    """
    from mlbspline.eval import evalMultivarSpline

    P = np.asarray(P, float).ravel(); T = np.asarray(T, float).ravel()
    grid = [P, T] + ([np.asarray(m, float).ravel()] if m is not None else [])
    nd = len(grid)
    x = np.empty(nd, dtype=object)
    for i, g in enumerate(grid):
        x[i] = g
    mesh = [g.reshape([-1 if i == j else 1 for j in range(nd)]) for i, g in enumerate(grid)]
    shape = tuple(g.size for g in grid)
    f = lambda *a: np.broadcast_to(np.asarray(G(*a), float), shape)

    def der(*d):
        return evalMultivarSpline(sp, x, list(d) + [0] * (nd - len(d)), allowExtrapolations=True)

    def shifted(i, h):
        a = list(mesh); a[i] = a[i] + h
        return f(*a)

    hP, hT = 1e-2, 1e-2
    G0 = f(*mesh)
    dGdP = (shifted(0, hP) - shifted(0, -hP)) / (2 * hP)
    d2GdT2 = (shifted(1, hT) - 2 * G0 + shifted(1, -hT)) / hT**2
    Tm = mesh[1]
    err = {
        'G': float(np.max(np.abs(der(0, 0) - G0))),
        'rho': float(np.max(np.abs(1e6 / der(1, 0) - 1e6 / dGdP))),
        'Cp': float(np.max(np.abs(-Tm * der(0, 2) - (-Tm * d2GdT2)))),
    }
    if nd == 3 and MW is not None:
        MWv = 18.01528e-3
        mm = mesh[2]
        ff = 1.0 + MW * mm
        hm = np.minimum(1e-4, 0.25 * np.where(mm > 0, mm, 1e-4))
        up = list(mesh); up[2] = mm + hm
        dn = list(mesh); dn[2] = np.maximum(mm - hm, 0.0)
        dGdm = (f(*up) - f(*dn)) / (up[2] - dn[2])
        muw_true = MWv * (G0 - ff * mm * dGdm)
        muw_spl = MWv * (der(0, 0, 0) - ff * mm * der(0, 0, 1))
        err['muw'] = float(np.max(np.abs(muw_spl - muw_true)))
    return err
