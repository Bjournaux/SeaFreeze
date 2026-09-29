"""
Thermodynamic properties from a Helmholtz-energy representation F(rho, T).

Counterpart of ``evalGibbs`` for equations of state stored as Helmholtz
splines; mirrors the Matlab ``fnFval.m`` / ``psi_val.m`` pair.

Two representations are supported, selected by ``sp['eos']``:

``'F_rhoT'``
    A plain tensor B-spline of F (J/kg) with knots[0] = rho (kg/m^3) and
    knots[1] = T (K).  Optional ``Tc``: the spline is in tau = log(T/Tc).

``'psi'``
    The psi-spline surface of lbf-thermo (JMB 2026): the residual
    dimensionless Helmholtz energy as a tensor B-spline in stretched
    coordinates plus reference terms,

        phi(delta, tau) = phi0 + dphi_ref + phi_crit + phi_2s + delta * psi(x, y)
        x = ln(rho/rhoc)/3,  y = ln(T/Tc),  delta = rho/rhoc,  tau = Tc/T
        F = R T phi

    with the Planck-Einstein ideal gas phi0, a reacting-mixture reference
    table (clamped at its edges, low-density extension joined C^1 below its
    floor), the tabulated KW2000 critical term (zero outside its box), an
    analytic low-T two-structure term, and a virial continuation of psi
    below its lowest density knot.  This is a field-by-field port of
    psiH2O_val.m.

Thermodynamics (SI inside; P and moduli returned in MPa)::

    P     = rho^2 F_r               G  = F + rho F_r       A = F
    S     = -F_T                    U  = F + T S           H = U + P/rho
    Cv    = -T F_TT                 Kt = rho dP/drho
    alpha = (dP/dT)_rho / Kt        Cp = Cv + T (dP/dT)^2 / (rho^2 dP/drho)
    Ks    = Kt Cp/Cv                vel = sqrt(Ks/rho)
    Kp    = 1 + rho (d2P/drho2) / (dP/drho)
    dP/drho = 2 rho F_r + rho^2 F_rr,  dP/dT = rho^2 F_rT,
    d2P/drho2 = 2 F_r + 4 rho F_rr + rho^2 F_rrr

Public API
----------
evalHelmholtzGrid(sp, PTm, *props, rhoT=False, allowExtrapolations=False)
evalHelmholtzScatter(sp, PTm, *props, rhoT=False, allowExtrapolations=False)
ideal_gas_props(sp, X, T, *props, rhoT=False)   ideal-gas part alone (psi)

``PTm`` follows the seafreeze conventions: a grid is ``np.array([X_vec, T_vec],
dtype=object)``, a scatter set is a 1-D object array of ``(X, T)`` tuples.
``X`` is pressure (MPa) by default, or density (kg/m^3) when ``rhoT=True``.
Output attributes use the lbftd names (G, S, U, H, A, rho, Cp, Cv, Kt, Kp,
Ks, alpha, vel) plus P (MPa) and T (K); the input coordinate is echoed.

Baptiste Journaux - 2026
"""

import json
import re
from types import SimpleNamespace

import numpy as np
from scipy.interpolate import NdBSpline

# lbftd property names supported here (pure phases only)
SUPPORTED = ('G', 'S', 'U', 'H', 'A', 'rho', 'Cp', 'Cv', 'Kt', 'Kp', 'Ks',
             'alpha', 'vel', 'P', 'T')
# ideal_gas_props names: getProp's water3 output set
IDEAL_GAS_PROPS = SUPPORTED + ('Js', 'gamma_Gruneisen')
_DERIV_NAMES = ('F', 'Fr', 'Frr', 'Frrr', 'FT', 'FTT', 'FrT')

# IAPWS-95 R6-95(2018) ideal-gas coefficients (fallback when the surface
# carries none)
_N0_IAPWS = np.array([-8.3204464837497, 6.6832105275932, 3.00632, 0.012436,
                      0.97315, 1.27950, 0.96956, 0.24873])
_G0_IAPWS = np.array([1.28728967, 3.53734222, 7.74073708, 9.24437796, 27.5075105])


# ---------------------------------------------------------------------------
# Spline-dict helpers
# ---------------------------------------------------------------------------
def _unwrap(v):
    """Recursively convert a Matlab struct as loaded by hdf5storage/scipy
    (structured ndarray) into a plain dict of ndarrays / floats / strings."""
    if isinstance(v, dict):
        return {k: _unwrap(x) for k, x in v.items()}
    if isinstance(v, np.void) and v.dtype.names:
        return {n: _unwrap(v[n]) for n in v.dtype.names}
    if isinstance(v, np.ndarray) and v.dtype.names:
        if v.shape == (1, 1) or v.size == 1:
            return {n: _unwrap(v[n].ravel()[0]) for n in v.dtype.names}
        return [_unwrap(e) for e in v.ravel()]
    if isinstance(v, np.ndarray) and v.dtype == object:
        return [_unwrap(e) for e in v.ravel()]
    if isinstance(v, np.ndarray) and v.dtype.kind in 'US':
        return str(v.ravel()[0]) if v.size == 1 else [str(e) for e in v.ravel()]
    if isinstance(v, np.ndarray) and v.size == 1 and v.dtype.kind in 'fiub':
        return float(v.ravel()[0])
    if isinstance(v, np.ndarray) and v.dtype.kind in 'US':
        return str(v)
    return v


def _field(sp, name):
    """A (sub)struct or scalar field of the spline dict, unwrapped once and cached."""
    cache = sp.setdefault('_helm_cache', {})
    key = 'field:' + name
    if key not in cache:
        cache[key] = _unwrap(sp[name]) if name in sp else None
    return cache[key]


def _scalar(sp, name, default=None):
    v = sp.get(name, default)
    if v is None:
        return default
    v = _unwrap(v)
    if isinstance(v, list):
        v = v[0]
    return float(v)


def _string(sp, name):
    v = _unwrap(sp.get(name))
    if v is None or isinstance(v, (float, list, dict)):
        return None
    return str(v)


def _knots(spd):
    """Knot vectors of a spline dict as 1-D float arrays."""
    ks = spd['knots']
    if isinstance(ks, np.ndarray) and ks.dtype == object:
        ks = list(ks.ravel())
    return [np.asarray(k, float).ravel() for k in ks]


def _orders(spd):
    o = spd['order']
    if isinstance(o, list):
        o = np.array(o)
    return [int(v) for v in np.asarray(o, float).ravel()]


def _ndspline(spd, cache_holder=None, tag=''):
    """scipy NdBSpline for a 2-D B-form spline dict (extrapolate=False: NaN
    outside the knot span, turned into 0 by callers to mirror Matlab fnval)."""
    holder = cache_holder if cache_holder is not None else spd
    cache = holder.setdefault('_helm_cache', {}) if isinstance(holder, dict) else {}
    key = 'nd:' + tag
    if key in cache:
        return cache[key]
    tx, ty = _knots(spd)
    kx, ky = _orders(spd)
    c = np.asarray(spd['coefs'], float)
    c = c.reshape(len(tx) - kx, len(ty) - ky)
    nd = NdBSpline((tx, ty), c, k=(kx - 1, ky - 1), extrapolate=False)
    cache[key] = nd
    return nd


def _ev(nd, pts, i, j):
    """Spline derivative (i, j) at points (n, 2); 0 outside the span."""
    if pts.shape[0] == 0:
        return np.zeros(0)
    return np.nan_to_num(nd(pts, nu=(i, j)))


def is_psi(sp):
    return _string(sp, 'eos') == 'psi'


def is_helmholtz(sp):
    return _string(sp, 'eos') in ('psi', 'F_rhoT')


def domain(sp):
    """((rho_lo, rho_hi), (T_lo, T_hi)) covered by the spline."""
    tx, ty = _knots(sp)
    if is_psi(sp):
        rhoc, Tc = _scalar(sp, 'rhoc'), _scalar(sp, 'Tc')
        # below the lowest density knot the psi surface is continued (virial)
        return (0.0, rhoc * np.exp(3 * tx[-1])), (Tc * np.exp(ty[0]), Tc * np.exp(ty[-1]))
    Tc = _scalar(sp, 'Tc')
    if Tc is not None:
        return (tx[0], tx[-1]), (Tc * np.exp(ty[0]), Tc * np.exp(ty[-1]))
    return (tx[0], tx[-1]), (ty[0], ty[-1])


# ---------------------------------------------------------------------------
# F derivatives: plain F(rho,T) B-spline
# ---------------------------------------------------------------------------
def _spline_derivs(sp, rho, T, need):
    Tc = _scalar(sp, 'Tc')
    t = np.log(T / Tc) if Tc is not None else T
    pts = np.stack([rho, t], -1)
    nd = _ndspline(sp, tag='main')
    d = {}
    if need['F']:
        d['F'] = nd(pts, nu=(0, 0))
    if need['Fr']:
        d['Fr'] = nd(pts, nu=(1, 0))
    if need['Frr']:
        d['Frr'] = nd(pts, nu=(2, 0))
    if need['Frrr']:
        d['Frrr'] = nd(pts, nu=(3, 0))
    need_FT = need['FT'] or (Tc is not None and need['FTT'])
    if need_FT:
        d['FT'] = nd(pts, nu=(0, 1))
    if need['FTT']:
        d['FTT'] = nd(pts, nu=(0, 2))
    if need['FrT']:
        d['FrT'] = nd(pts, nu=(1, 1))
    if Tc is not None:
        # d/dT = (1/T) d/dtau ; d2/dT2 = (d2/dtau2 - d/dtau) / T^2
        if need['FTT']:
            d['FTT'] = (d['FTT'] - d['FT']) / T ** 2
        if need_FT:
            d['FT'] = d['FT'] / T
        if need['FrT']:
            d['FrT'] = d['FrT'] / T
    if need_FT and not need['FT']:
        del d['FT']
    return d


# ---------------------------------------------------------------------------
# F derivatives: psi surface (port of psi_val.m / psiH2O_val.m)
# ---------------------------------------------------------------------------
def _lowT2s_params(txt):
    p = {}
    for k, v in re.findall(r'"(\w+)"\s*:\s*([^,}]+)', txt):
        v = v.strip()
        p[k] = None if v == 'null' else float(v)
    return p


def _lowT2s_phi(d, tau, p, R):
    """phi_2s = min_x [x G + x ln x + (1-x) ln(1-x) + omega x (1-x)] and its
    delta/tau derivatives at the equilibrium x (twostructure_lowT.TwoStructureLowT)."""
    TC = 647.096
    RHOC = 322.0
    if p.get('window') is not None:
        raise ValueError('windowed low-T two-structure term not supported')
    s_ = p['s']
    th = p['Tstar'] / TC
    be = p['b'] / (R * TC)
    ds = p['rhostar'] / RHOC
    la = p.get('lam') or 0.0
    L = np.log(d / ds)
    if p.get('rho_m') is not None:
        Lm = np.log(p['rho_m'] / p['rhostar'])
        Q = L - L ** 2 / (2 * Lm)
        Q1 = 1 - L / Lm
        Q2 = -1 / Lm + 0 * L
    else:
        Q = L
        Q1 = 1 + 0 * L
        Q2 = 0 * L
    g = s_ * (th * tau - 1) + la * (1 / (th * tau) - 1) + be * tau * Q
    gd = be * tau * Q1 / d
    gdd = be * tau * (Q2 - Q1) / d ** 2
    gt = s_ * th - la / (th * tau ** 2) + be * Q
    gtt = 2 * la / (th * tau ** 3)
    gdt = be * Q1 / d
    Lw = np.log(d / (p['rho_top'] / RHOC))
    om_h = p.get('omega_h')
    if om_h is not None and om_h > 0:
        h = om_h
        zz = -Lw / h
        sp_ = h * (np.maximum(zz, 0) + np.log1p(np.exp(-np.abs(zz))))
        sg = 1 / (1 + np.exp(-zz))
        om = p['omega_top'] - p['omega1'] * sp_
        dL = p['omega1'] * sg
        d2L = -(p['omega1'] / h) * sg * (1 - sg)
        omd = dL / d
        omdd = (d2L - dL) / d ** 2
    else:
        om = p['omega_top'] + p['omega1'] * Lw
        omd = p['omega1'] / d
        omdd = -p['omega1'] / d ** 2
    # x from G + u + omega (1 - 2x) = 0, u = ln(x/(1-x)): bisection + Newton
    a = np.abs(om) + 1
    lo = np.clip(-g - a, -700, 700)
    hi = np.clip(-g + a, -700, 700)
    for _ in range(34):
        mid = 0.5 * (lo + hi)
        pos = (g + mid + om * (1 - 2 / (1 + np.exp(-mid)))) > 0
        hi = np.where(pos, mid, hi)
        lo = np.where(pos, lo, mid)
    u = 0.5 * (lo + hi)
    for _ in range(3):
        xx = 1 / (1 + np.exp(-u))
        fu = g + u + om * (1 - 2 * xx)
        u = np.clip(u - fu / (1 - 2 * om * xx * (1 - xx)), -700, 700)
    xx = 1 / (1 + np.exp(-u))
    qq = np.maximum(xx * (1 - xx), 1e-300)
    ent = -(np.maximum(u, 0) + np.log1p(np.exp(-np.abs(u)))) + xx * u
    Pxx = 1 / qq - 2 * om
    Pxd = gd + omd * (1 - 2 * xx)
    Pxt = gt
    return {'p': xx * g + ent + om * qq,
            'd': xx * gd + omd * qq,
            't': xx * gt,
            'dd': xx * gdd + omdd * qq - Pxd ** 2 / Pxx,
            'tt': xx * gtt - Pxt ** 2 / Pxx,
            'dt': xx * gdt - Pxd * Pxt / Pxx}


def _psi_derivs(sp, rho, T, need):
    n = rho.size
    Tc, rhoc, R = _scalar(sp, 'Tc'), _scalar(sp, 'rhoc'), _scalar(sp, 'R')
    ref = _field(sp, 'dphi_ref')
    if any(k in sp for k in ('twostate_params', 'twostate_s_params')) or \
       (ref is not None and any(k in ref for k in ('twostate_params', 'twostate_s_params'))):
        raise ValueError('retired two-state terms are not supported')
    if ref is not None and 'crit_params' in ref and 'crit_table' not in ref:
        raise ValueError('direct KW2000 critical term is not supported; the surface must carry crit_table')

    if need['Frrr']:
        # central finite difference of Frr (the reference terms carry no
        # analytic third derivative)
        h = 1e-4
        nd_ = dict.fromkeys(_DERIV_NAMES, False)
        nd_['Frr'] = True
        dp = _psi_derivs(sp, rho * (1 + h), T, nd_)
        dm = _psi_derivs(sp, rho * (1 - h), T, nd_)
        need = dict(need)
        need['Frrr'] = False
        d = _psi_derivs(sp, rho, T, need)
        d['Frrr'] = (dp['Frr'] - dm['Frr']) / (2 * h * rho)
        return d

    wp = need['F'] or need['FT']
    wd = need['Fr'] or need['FrT']
    wdd = need['Frr']
    wt = need['FT']
    wtt = need['FTT']
    wdt = need['FrT']
    wy = wt or wtt or wdt

    with np.errstate(divide='ignore', invalid='ignore'):
        x = np.log(rho / rhoc) / 3
        y = np.log(T / Tc)
    tx, ty = _knots(sp)
    ok = np.isfinite(x) & np.isfinite(y) & (rho > 0) & (x <= tx[-1]) & (y >= ty[0]) & (y <= ty[-1])

    delta = rho / rhoc
    tau = Tc / T
    z = np.full(n, np.nan)
    phir, phir_d, phir_dd, phir_t, phir_tt, phir_dt = (z.copy() for _ in range(6))

    if ok.any():
        xo, yo, do, to = x[ok], y[ok], delta[ok], tau[ok]
        pts = np.stack([xo, yo], -1)
        nd = _ndspline(sp, tag='main')

        # ---- spline residual psi and derivatives in (x, y) ----------------
        s = _ev(nd, pts, 0, 0)
        sx = _ev(nd, pts, 1, 0)
        sxx = _ev(nd, pts, 2, 0) if wdd else None
        sy = _ev(nd, pts, 0, 1) if wy else None
        syy = _ev(nd, pts, 0, 2) if wtt else None
        sxy = _ev(nd, pts, 1, 1) if wdt else None

        # virial continuation below the lowest density knot (C^1 join)
        x0 = tx[0]
        below = xo < x0
        if below.any():
            qb = np.stack([np.full(below.sum(), x0), yo[below]], -1)
            r = np.exp(3 * (xo[below] - x0))
            s00, s10 = _ev(nd, qb, 0, 0), _ev(nd, qb, 1, 0)
            s[below] = s00 + s10 * (r - 1) / 3
            sx[below] = s10 * r
            if wdd:
                sxx[below] = 3 * s10 * r
            if wy:
                s01, s11 = _ev(nd, qb, 0, 1), _ev(nd, qb, 1, 1)
                sy[below] = s01 + s11 * (r - 1) / 3
                if wtt:
                    s02, s12 = _ev(nd, qb, 0, 2), _ev(nd, qb, 1, 2)
                    syy[below] = s02 + s12 * (r - 1) / 3
                if wdt:
                    sxy[below] = s11 * r

        pr = do * s
        pr_d = s + sx / 3
        pr_dd = (sx / 3 + sxx / 9) / do if wdd else None
        pr_t = -do * sy / to if wt else None
        pr_tt = do * (syy + sy) / to ** 2 if wtt else None
        pr_dt = -(sy + sxy / 3) / to if wdt else None

        # ---- reference table dphi_ref (clamped; low-density extension) ---
        if ref is not None:
            rkx, rky = _knots(ref)
            ndr = _ndspline(ref, cache_holder=sp, tag='ref')
            Xc = np.clip(xo, rkx[0], rkx[-1])
            Yc = np.clip(yo, rky[0], rky[-1])
            oob = (xo != Xc) | (yo != Yc)
            pr2 = np.stack([Xc, Yc], -1)

            def z0(v):
                v[oob] = 0.0
                return v
            f = _ev(ndr, pr2, 0, 0)
            fx = z0(_ev(ndr, pr2, 1, 0))
            fxx = z0(_ev(ndr, pr2, 2, 0)) if wdd else None
            fy = z0(_ev(ndr, pr2, 0, 1)) if wy else None
            fyy = z0(_ev(ndr, pr2, 0, 2)) if wtt else None
            fxy = z0(_ev(ndr, pr2, 1, 1)) if wdt else None
            if 'low' in ref:
                low = ref['low']
                lkx, _ = _knots(low)
                ndl = _ndspline(low, cache_holder=sp, tag='low')
                bl = (xo < rkx[0]) & (xo >= lkx[0]) & (yo >= rky[0]) & (yo <= rky[-1])
                if bl.any():
                    ql = np.stack([xo[bl], yo[bl]], -1)
                    q0 = np.stack([np.full(bl.sum(), rkx[0]), yo[bl]], -1)
                    r = np.exp(3 * (xo[bl] - rkx[0]))

                    def Lv(i, j):
                        return _ev(ndl, ql, i, j)

                    def D(i, j):
                        return _ev(ndr, q0, i, j) - _ev(ndl, q0, i, j)
                    D10 = D(1, 0)
                    f[bl] = Lv(0, 0) + D(0, 0) + D10 * (r - 1) / 3
                    fx[bl] = Lv(1, 0) + D10 * r
                    if wdd:
                        fxx[bl] = Lv(2, 0) + 3 * D10 * r
                    if wy:
                        D01, D11 = D(0, 1), D(1, 1)
                        fy[bl] = Lv(0, 1) + D01 + D11 * (r - 1) / 3
                        if wtt:
                            fyy[bl] = Lv(0, 2) + D(0, 2) + D(1, 2) * (r - 1) / 3
                        if wdt:
                            fxy[bl] = Lv(1, 1) + D11 * r
            pr = pr + f
            pr_d = pr_d + fx / (3 * do)
            if wdd:
                pr_dd = pr_dd + (fxx / 9 - fx / 3) / do ** 2
            if wt:
                pr_t = pr_t - fy / to
            if wtt:
                pr_tt = pr_tt + (fyy + fy) / to ** 2
            if wdt:
                pr_dt = pr_dt - fxy / (3 * do * to)

        # ---- tabulated KW2000 critical term ---------------------------------
        ct = _field(sp, 'crit_table')
        if ct is None and ref is not None:
            ct = ref.get('crit_table')
        if ct is not None:
            c0 = float(ct['c0'])
            ckx, cky = _knots(ct)
            ndc = _ndspline(ct, cache_holder=sp, tag='crit')
            xq = 1 - to
            yq = do - 1
            inb = (xq >= ckx[0]) & (xq <= ckx[-1]) & (yq >= cky[0]) & (yq <= cky[-1])
            if inb.any():
                pq = np.stack([xq[inb], yq[inb]], -1)
                m = ok.sum()
                g, gy, gyy, gx, gxx, gxy = (np.zeros(m) for _ in range(6))
                g[inb] = _ev(ndc, pq, 0, 0)
                gy[inb] = _ev(ndc, pq, 0, 1)
                if wdd:
                    gyy[inb] = _ev(ndc, pq, 0, 2)
                if wy:
                    gx[inb] = _ev(ndc, pq, 1, 0)
                if wtt:
                    gxx[inb] = _ev(ndc, pq, 2, 0)
                if wdt:
                    gxy[inb] = _ev(ndc, pq, 1, 1)
                pr = pr + c0 * g / do
                pr_d = pr_d + c0 * (gy / do - g / do ** 2)
                if wdd:
                    pr_dd = pr_dd + c0 * (gyy / do - 2 * gy / do ** 2 + 2 * g / do ** 3)
                if wt:
                    pr_t = pr_t - c0 * gx / do
                if wtt:
                    pr_tt = pr_tt + c0 * gxx / do
                if wdt:
                    pr_dt = pr_dt - c0 * (gxy / do - gx / do ** 2)

        # ---- low-T two-structure term (analytic) ----------------------------
        lj = _string(sp, 'lowT2s_params')
        if lj is None and ref is not None and 'lowT2s_params' in ref:
            lj = str(ref['lowT2s_params'])
        if lj:
            q = _lowT2s_phi(do, to, _lowT2s_params(lj), R)
            pr = pr + q['p']
            pr_d = pr_d + q['d']
            if wdd:
                pr_dd = pr_dd + q['dd']
            if wt:
                pr_t = pr_t + q['t']
            if wtt:
                pr_tt = pr_tt + q['tt']
            if wdt:
                pr_dt = pr_dt + q['dt']

        phir[ok] = pr
        phir_d[ok] = pr_d
        if wdd:
            phir_dd[ok] = pr_dd
        if wt:
            phir_t[ok] = pr_t
        if wtt:
            phir_tt[ok] = pr_tt
        if wdt:
            phir_dt[ok] = pr_dt

    # ---- ideal-gas part ------------------------------------------------------
    phi0, phi0_t, phi0_tt = _phi0(sp, delta, tau)
    with np.errstate(divide='ignore', invalid='ignore'):
        # ---- F and its (rho,T) derivatives ------------------------------------
        #   F = R T phi;  d/drho = (1/rhoc) d/ddelta;  d tau/dT = -tau/T
        RT = R * T
        d = {}
        if wp:
            phi = phi0 + phir
            if need['F']:
                d['F'] = RT * phi
            if need['FT']:
                d['FT'] = R * (phi - tau * (phi0_t + phir_t))
        if wd:
            phi_d = 1 / delta + phir_d
            if need['Fr']:
                d['Fr'] = RT * phi_d / rhoc
            if need['FrT']:
                d['FrT'] = R / rhoc * (phi_d - tau * phir_dt)
        if need['Frr']:
            d['Frr'] = RT * (-1 / delta ** 2 + phir_dd) / rhoc ** 2
        if need['FTT']:
            d['FTT'] = R * tau ** 2 * (phi0_tt + phir_tt) / T
    for k in d:
        d[k] = np.where(ok, d[k], np.nan)
    return d


def _phi0(sp, delta, tau):
    """Planck-Einstein ideal-gas part phi0 and its tau derivatives."""
    if 'phi0_n0' in sp and 'phi0_g0' in sp:
        n0 = np.asarray(_unwrap(sp['phi0_n0']), float).ravel()
        g0 = np.asarray(_unwrap(sp['phi0_g0']), float).ravel()
    else:
        n0, g0 = _N0_IAPWS, _G0_IAPWS
    with np.errstate(divide='ignore', invalid='ignore'):
        phi0 = np.log(delta) + n0[0] + n0[1] * tau + n0[2] * np.log(tau)
        phi0_t = n0[1] + n0[2] / tau
        phi0_tt = -n0[2] / tau ** 2
        for i, gi in enumerate(g0):
            e = np.exp(-gi * tau)
            phi0 = phi0 + n0[3 + i] * np.log(1 - e)
            phi0_t = phi0_t + n0[3 + i] * gi * (1 / (1 - e) - 1)
            phi0_tt = phi0_tt - n0[3 + i] * gi ** 2 * e / (1 - e) ** 2
    return phi0, phi0_t, phi0_tt


def ideal_gas_props(sp, X, T, *props, rhoT=False):
    """Properties of the surface's ideal-gas part alone (Z = 1).  psi surfaces only.

    This is the dilute-vapour limit of the EOS; it needs no spline and so
    stays defined below the surface's lowest temperature knot.  Used by
    seafreeze.sublimation(..., dilute_extension=True) and phase_map.  With
    the Planck-Einstein phi0(delta, tau) and its tau derivatives::

        rho = P/(R T)                 A = R T phi0          G = A + R T
        S = R (tau phi0_t - phi0)     U = R T tau phi0_t    H = U + R T
        Cv = -R tau^2 phi0_tt         Cp = Cv + R           alpha = 1/T
        Kt = P    Ks = Kt Cp/Cv    Kp = 1    vel = sqrt(Cp/Cv R T)
        Js = T alpha/(rho Cp) = 1/(rho Cp) (x 1e6, K/MPa)    gamma_Gruneisen = R/Cv

    Js is getProp's isentropic dT/dP, not the Joule-Thomson coefficient
    (which vanishes for an ideal gas).  Above ~700 K the surface's reference
    term (reacting mixture) makes the dilute fluid depart from this ideal gas.

    :param X:     P (MPa), or rho (kg/m^3) with rhoT=True; broadcasts with T
    :param T:     temperature (K)
    :param props: property names, as getProp's water3 output (same units);
                  none = all
    :return:      SimpleNamespace of arrays with the broadcast shape of X, T
    """
    if not is_psi(sp):
        raise ValueError('ideal_gas_props needs a psi surface (it carries phi0)')
    if not props:
        props = IDEAL_GAS_PROPS
    unknown = set(props) - set(IDEAL_GAS_PROPS)
    if unknown:
        raise ValueError('unsupported property name(s): ' + ', '.join(sorted(unknown)))
    X, T = (np.array(a, float) for a in np.broadcast_arrays(np.asarray(X, float), np.asarray(T, float)))
    R, rhoc, Tc = _scalar(sp, 'R'), _scalar(sp, 'rhoc'), _scalar(sp, 'Tc')
    RT = R * T
    with np.errstate(divide='ignore', invalid='ignore'):
        if rhoT:
            rho = X
            P = rho * RT / 1e6
        else:
            P = X
            rho = P * 1e6 / RT
        tau = Tc / T
        phi0, phi0_t, phi0_tt = _phi0(sp, rho / rhoc, tau)
        U = RT * tau * phi0_t
        Cv = -R * tau ** 2 * phi0_tt
        Cp = Cv + R
        alpha = 1 / T
        out = {'G': RT * (phi0 + 1.0), 'S': R * (tau * phi0_t - phi0), 'U': U, 'H': U + RT,
               'A': RT * phi0, 'rho': rho, 'Cp': Cp, 'Cv': Cv, 'Kt': P,
               'Kp': np.where(np.isnan(rho), np.nan, 1.0), 'Ks': P * Cp / Cv,
               'alpha': alpha, 'vel': np.sqrt(Cp / Cv * RT),
               'Js': T * alpha / (rho * Cp) * 1e6,                   # K/MPa
               'gamma_Gruneisen': R / Cv, 'P': P, 'T': T}
    return SimpleNamespace(**{k: out[k] for k in props})


def ideal_gas(sp, P, T):
    """Gibbs energy (J/kg) and density (kg/m^3) of the surface's ideal-gas
    part alone (Z = 1) at P (MPa), T (K).  psi surfaces only; see
    ideal_gas_props for the full property set.
    """
    o = ideal_gas_props(sp, P, T, 'G', 'rho')
    return o.G, o.rho


# ---------------------------------------------------------------------------
# Dispatch, P(rho,T) and its inversion
# ---------------------------------------------------------------------------
def _helm_derivs(sp, rho, T, need):
    rho = np.asarray(rho, float).ravel()
    T = np.asarray(T, float).ravel()
    nd_ = dict.fromkeys(_DERIV_NAMES, False)
    nd_.update(need)
    if rho.size == 0:
        return {k: np.zeros(0) for k, v in nd_.items() if v}
    if is_psi(sp):
        return _psi_derivs(sp, rho, T, nd_)
    return _spline_derivs(sp, rho, T, nd_)


def _P_dPdrho(sp, rho, T):
    """P (Pa) and (dP/drho)_T (Pa m^3/kg) at scattered points."""
    rho = np.asarray(rho, float).ravel()
    d = _helm_derivs(sp, rho, T, {'Fr': True, 'Frr': True})
    return rho ** 2 * d['Fr'], 2 * rho * d['Fr'] + rho ** 2 * d['Frr']


def _newton_bracket(sp, a, b, Pt, ta):
    """Bracket-safeguarded Newton for rho^2 F_r = Pt (Pa), P increasing on [a, b]."""
    a = a.copy(); b = b.copy()
    pa, _ = _P_dPdrho(sp, a, ta)
    pb, _ = _P_dPdrho(sp, b, ta)
    with np.errstate(divide='ignore', invalid='ignore'):
        x = a + (Pt - pa) * (b - a) / (pb - pa)
    bad = ~np.isfinite(x) | (x < a) | (x > b)
    x[bad] = 0.5 * (a[bad] + b[bad])
    live = np.ones(x.size, dtype=bool)
    for _ in range(80):
        f, df = _P_dPdrho(sp, x[live], ta[live])
        f = f - Pt[live]
        al, bl, xl = a[live], b[live], x[live]
        al = np.where(f < 0, xl, al)
        bl = np.where(f > 0, xl, bl)
        with np.errstate(divide='ignore', invalid='ignore'):
            xn = xl - f / df
        out = ~np.isfinite(xn) | (xn <= al) | (xn >= bl)
        xn = np.where(out, 0.5 * (al + bl), xn)
        step = np.abs(xn - xl)
        a[live], b[live], x[live] = al, bl, xn
        done = (step <= 1e-12 * xn) | (f == 0)
        live[np.flatnonzero(live)[done]] = False
        if not live.any():
            break
    return x


def _invert_P(sp, P, T, rlim, Tlim, branch='stable'):
    """Density solving rho^2 dF/drho = P (P in MPa) at each point.

    Brackets every root on a coarse density grid (per distinct temperature)
    where P increases with rho (mechanically stable), refines the least and
    the most dense of them, and returns the one selected by ``branch``:
    'stable' (lower Gibbs energy), 'liquid' (densest) or 'vapor' (least dense).
    NaN where no stable root exists.
    """
    if branch not in ('stable', 'liquid', 'vapor'):
        raise ValueError("branch must be 'stable', 'liquid' or 'vapor'")
    P = np.asarray(P, float).ravel()
    T = np.asarray(T, float).ravel()
    n = P.size
    rho = np.full(n, np.nan)
    inT = np.isfinite(P) & np.isfinite(T) & (T >= Tlim[0]) & (T <= Tlim[1])
    if not inT.any():
        return rho
    Ppa = P * 1e6

    nr = 600
    if is_psi(sp):
        # geometric spacing, extended below the lowest density knot into the
        # virial continuation (dilute vapour)
        tx, _ = _knots(sp)
        rhoc = _scalar(sp, 'rhoc')
        x_lo = min(tx[0], np.log(1e-12 / rhoc) / 3)
        rg = rhoc * np.exp(3 * np.linspace(x_lo, tx[-1], nr))
        # plus a fine uniform sampling of the dense fluid, so that the small
        # (dP/drho)_T loops inside the dome near Tc never share a bracket
        # with the true liquid or vapour root
        rtop = rhoc * np.exp(3 * tx[-1])
        rg = np.unique(np.r_[rg, np.arange(1.0, min(2000.0, rtop), 2.0)])
        nr = rg.size
    else:
        rg = np.linspace(rlim[0], rlim[1], nr)
    klo = np.zeros(n, dtype=int)
    khi = np.zeros(n, dtype=int)
    idx = np.flatnonzero(inT)
    Tu, jT = np.unique(T[idx], return_inverse=True)
    nTu = Tu.size
    pg, _ = _P_dPdrho(sp, np.tile(rg, nTu), np.repeat(Tu, nr))
    pg = pg.reshape(nTu, nr)
    kidx = np.arange(1, nr)[:, None]
    for j in range(nTu):
        pts = idx[jT == j]
        col = pg[j]
        s = col[:, None] - Ppa[pts][None, :]                    # nr x npts
        up = (s[:-1] <= 0) & (s[1:] >= 0) & (np.diff(col) > 0)[:, None]
        khi[pts] = np.max(np.where(up, kidx, 0), axis=0)        # densest crossing
        kl = np.min(np.where(up, kidx, nr + 1), axis=0)          # least dense crossing
        klo[pts] = np.where(kl > nr, 0, kl)
    if branch == 'liquid':
        klo = khi.copy()
    elif branch == 'vapor':
        khi = klo.copy()

    ih = np.flatnonzero(khi > 0)
    if ih.size == 0:
        return rho
    rho[ih] = _newton_bracket(sp, rg[khi[ih] - 1], rg[khi[ih]], Ppa[ih], T[ih])
    two = np.flatnonzero((khi > 0) & (klo != khi))
    if two.size:
        rlo = _newton_bracket(sp, rg[klo[two] - 1], rg[klo[two]], Ppa[two], T[two])
        nd_ = {'F': True, 'Fr': True}
        dh = _helm_derivs(sp, rho[two], T[two], nd_)
        dl = _helm_derivs(sp, rlo, T[two], nd_)
        Gh = dh['F'] + rho[two] * dh['Fr']
        Gl = dl['F'] + rlo * dl['Fr']
        usev = (Gl < Gh) | ~np.isfinite(Gh)
        rho[two[usev]] = rlo[usev]
    return rho


# ---------------------------------------------------------------------------
# Property assembly
# ---------------------------------------------------------------------------
def _eval(sp, X, T, props, rhoT, gridded, branch='stable'):
    """Core evaluation.  X, T are 1-D (grid vectors) or equal-length (scatter)."""
    if not is_helmholtz(sp):
        raise ValueError("spline is not a Helmholtz representation (sp['eos'] must be 'psi' or 'F_rhoT')")
    if not props:
        props = tuple(SUPPORTED)
    unknown = set(props) - set(SUPPORTED)
    if unknown:
        raise ValueError('unsupported property name(s): ' + ', '.join(sorted(unknown)))
    want = {p: (p in props) for p in SUPPORTED}

    X = np.asarray(X, float)
    T = np.asarray(T, float)
    if gridded:
        Xm, Tm = np.meshgrid(X, T, indexing='ij')
    else:
        Xm, Tm = X, T
    sz = Xm.shape
    rlim, Tlim = domain(sp)

    need_Cv = want['Cv'] or want['Cp'] or want['Ks'] or want['vel']
    need_alpha = want['alpha'] or want['Cp'] or want['Ks'] or want['vel']
    need = {'F': want['G'] or want['U'] or want['H'] or want['A'],
            'Fr': True,
            'Frr': need_alpha or want['Kt'] or want['Kp'],
            'Frrr': want['Kp'],
            'FT': want['S'] or want['U'] or want['H'],
            'FTT': need_Cv,
            'FrT': need_alpha}

    if rhoT:
        rhom = Xm.astype(float).copy()
        rhom[(rhom < rlim[0]) | (rhom > rlim[1])] = np.nan
    else:
        rhom = _invert_P(sp, Xm.ravel(), Tm.ravel(), rlim, Tlim, branch).reshape(sz)
    rhom = np.where((Tm < Tlim[0]) | (Tm > Tlim[1]), np.nan, rhom)

    ok = np.isfinite(rhom)
    dd = _helm_derivs(sp, rhom[ok], Tm[ok], need)
    d = {}
    for k, v in dd.items():
        full = np.full(sz, np.nan)
        full[ok] = v
        d[k] = full

    r = rhom
    with np.errstate(divide='ignore', invalid='ignore'):
        Pp = r ** 2 * d['Fr']                                        # Pa
        if need['Frr']:
            dPdr = 2 * r * d['Fr'] + r ** 2 * d['Frr']               # Pa m^3/kg
        if need_Cv:
            Cv = -Tm * d['FTT']
        if need_alpha:
            dPdT = r ** 2 * d['FrT']                                 # Pa/K
            alpha = dPdT / (r * dPdr)
            if need_Cv:
                Cp = Cv + Tm * dPdT ** 2 / (r ** 2 * dPdr)
        if need['Frr']:
            Kt_Pa = r * dPdr

        out = {}
        if want['G']:
            out['G'] = d['F'] + Pp / r
        if want['S']:
            out['S'] = -d['FT']
        if want['U']:
            out['U'] = d['F'] - Tm * d['FT']
        if want['H']:
            out['H'] = d['F'] - Tm * d['FT'] + Pp / r
        if want['A']:
            out['A'] = d['F']
        if want['rho']:
            out['rho'] = X.copy() if (rhoT and gridded) else rhom
        if want['Cp']:
            out['Cp'] = Cp
        if want['Cv']:
            out['Cv'] = Cv
        if want['Kt']:
            out['Kt'] = Kt_Pa / 1e6
        if want['Ks']:
            out['Ks'] = Kt_Pa * Cp / Cv / 1e6
        if want['Kp']:
            d2Pdr2 = 2 * d['Fr'] + 4 * r * d['Frr'] + r ** 2 * d['Frrr']
            out['Kp'] = 1 + r * d2Pdr2 / dPdr
        if want['alpha']:
            out['alpha'] = alpha
        if want['vel']:
            w2 = Kt_Pa * Cp / Cv / r
            out['vel'] = np.sqrt(np.where(w2 < 0, np.nan, w2))       # NaN where unstable
        if want['P']:
            out['P'] = (X.copy() if gridded else Xm) if not rhoT else Pp / 1e6
        if want['T']:
            out['T'] = T.copy() if gridded else Tm
    return SimpleNamespace(**out)


def evalHelmholtzGrid(sp, PTm, *props, rhoT=False, branch='stable', allowExtrapolations=False):
    """Properties on the tensor grid PTm = array([X_vec, T_vec], dtype=object).

    X is P (MPa), or rho (kg/m^3) when rhoT=True.  Output arrays are
    (len(X), len(T)).  ``branch`` (P input only) picks the fluid root when
    both a vapour-like and a liquid-like stable root exist: 'stable' (lower
    Gibbs energy, default), 'liquid' or 'vapor' (metastable allowed).  Points outside the spline domain are NaN
    (allowExtrapolations is accepted for signature parity and ignored: a
    Helmholtz surface has no meaningful extrapolation).
    """
    return _eval(sp, np.asarray(PTm[0], float).ravel(), np.asarray(PTm[1], float).ravel(),
                 props, rhoT, True, branch)


def evalHelmholtzScatter(sp, PTm, *props, rhoT=False, branch='stable', allowExtrapolations=False):
    """Properties at scattered points PTm = object array of (X, T) tuples."""
    X = np.array([float(t[0]) for t in PTm])
    T = np.array([float(t[1]) for t in PTm])
    return _eval(sp, X, T, props, rhoT, False, branch)
