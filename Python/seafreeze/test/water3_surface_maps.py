"""Property maps of the Helmholtz fluid 'water3' over its whole (P,T) range.

    python3 -m seafreeze.test.water3_surface_maps OUTDIR

Writes water3_surface_maps.png: density, Cp, sound speed, thermal expansivity,
compressibility factor Z = P/(rho R T) and the Grueneisen parameter on a log P -
log T grid (1e-6 MPa - 10 TPa, 230 K - 150 kK), stable branch (vapour below the
saturation curve).  Overlays: water3's own saturation curve and critical point,
and the melting curve of the stable solid (psiEOS 'dq2026' model: IAPWS R14-08,
Datchi 2000, Queyroux 2020, French & Hamel) — below it the surface still returns
numbers but the stable phase is a solid (hatched).
"""
import os
import sys
import time
import warnings

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm, Normalize

import seafreeze as sf
from seafreeze.seafreeze import defpath

INK, INK2 = '#0b0b0b', '#52514e'
plt.rcParams.update({'font.size': 9.5, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'figure.facecolor': '#fcfcfb',
                     'axes.facecolor': '#ecebe7'})


# ---- melting curve of the stable solid (port of psiEOS.m melt_T, JMB 2026) ----
def _p_melt_Ih(T):
    th = T / 273.16
    return 611.657e-6 * (1 + 0.119539337e7 * (1 - th ** 3) + 0.808183159e5 * (1 - th ** 25.75)
                         + 0.333826860e4 * (1 - th ** 103.75))


def _p_melt_hp(T):
    T = np.asarray(T, float); p = np.full(T.shape, np.nan)
    m = (T >= 251.165) & (T < 256.164); th = T[m] / 251.165; p[m] = 208.566 * (1 - 0.299948 * (1 - th ** 60))
    m = (T >= 256.164) & (T < 273.31); th = T[m] / 256.164; p[m] = 350.1 * (1 - 1.18721 * (1 - th ** 8))
    m = (T >= 273.31) & (T < 355.0); th = T[m] / 273.31; p[m] = 632.4 * (1 - 1.07476 * (1 - th ** 4.6))
    m = (T >= 355.0) & (T <= 715.0); th = T[m] / 355.0
    p[m] = 2216.0 * np.exp(0.173683e1 * (1 - 1 / th) - 0.544606e-1 * (1 - th ** 5) + 0.806106e-7 * (1 - th ** 22))
    t7 = 715 / 355
    P715 = 2216.0 * np.exp(0.173683e1 * (1 - 1 / t7) - 0.544606e-1 * (1 - t7 ** 5) + 0.806106e-7 * (1 - t7 ** 22))
    n = np.log(45000 / P715) / np.log(1600 / 715)
    m = T > 715.0; p[m] = P715 * (T[m] / 715) ** n
    return p


def melt_T(P):
    """Melting temperature (K) of the stable solid at P (MPa); psiEOS.m melt_T."""
    P = np.asarray(P, float).ravel(); Tm = np.full(P.shape, np.nan)
    PT_si = np.array([[45, 850 * ((45 - 14.6) / 3.44 + 1) ** (1 / 4.33)], [330.6, 6000], [627.6, 7000],
                      [1034.2, 8000], [1592.9, 10000], [3734.6, 10000], [4664.8, 12000], [9311.3, 12000]])
    Tm[P < 611.657e-6] = 273.16
    ih = (P >= 611.657e-6) & (P < 208.566)
    a = np.full(ih.sum(), 251.165); b = np.full(ih.sum(), 273.16); Pi = P[ih]
    for _ in range(60):
        m = 0.5 * (a + b); up = _p_melt_Ih(m) > Pi; a = np.where(up, m, a); b = np.where(up, b, m)
    Tm[ih] = 0.5 * (a + b)
    hp = (P >= 208.566) & (P < 2170)
    a = np.full(hp.sum(), 251.165); b = np.full(hp.sum(), 2e4); Ph = P[hp]
    for _ in range(80):
        m = 0.5 * (a + b); up = _p_melt_hp(m) < Ph; a = np.where(up, m, a); b = np.where(up, b, m)
    Tm[hp] = 0.5 * (a + b)
    vii = (P >= 2170) & (P <= 45000)
    Pg = P[vii] / 1e3
    Td = 354.8 * np.maximum((Pg - 2.17) / 1.253 + 1, 1e-12) ** (1 / 3)
    xq = (Pg - 14.6) / 3.44 + 1
    Tq = np.where(xq > 0, 850 * np.maximum(xq, 0) ** (1 / 4.33), 0)
    Tm[vii] = np.maximum(Td, Tq)
    si = P > 45000
    Tm[si] = np.interp(np.log(np.minimum(P[si] / 1e3, PT_si[-1, 0])), np.log(PT_si[:, 0]), PT_si[:, 1])
    return Tm


def make(outdir):
    P = np.geomspace(1e-6, 1e7, 260)
    T = np.geomspace(230, 1.5e5, 220)
    t0 = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        o = sf.getProp(np.array([P, T], dtype=object), 'water3', defpath,
                       'rho', 'Cp', 'vel', 'alpha', 'gamma_Gruneisen')
    print(f'evaluated {P.size * T.size} states in {time.time() - t0:.1f} s')
    R = 461.5231157260608
    Z = (P[:, None] * 1e6) / (o.rho * R * T[None, :])
    Tsat = np.linspace(230, 646.9, 300)
    sat = sf.saturation(Tsat)
    Pm = np.geomspace(1e-6, 1e7, 600)
    Tm = melt_T(Pm)
    # below the triple-point pressure the solid/vapour boundary is sublimation
    Tsub = np.linspace(230, 273.16, 60)
    psub = sf.sublimation(Tsub).P
    low = Pm < 611.657e-6
    Tm[low] = np.interp(np.log(Pm[low]), np.log(psub), Tsub, left=np.nan)

    panels = [
        ('Density ρ (kg/m³)', o.rho, LogNorm(1e-4, 1.6e4), 'Blues'),
        ('Isobaric heat capacity Cp (J/kg/K)', o.Cp, LogNorm(1.5e3, 1e5), 'Blues'),
        ('Sound speed (m/s)', o.vel, LogNorm(3e2, 3e4), 'Blues'),
        ('Thermal expansivity α (1/K)', o.alpha, SymLogNorm(1e-5, vmin=-3e-3, vmax=3e-3), 'RdBu_r'),
        ('Compressibility factor Z = P/(ρRT)', Z, LogNorm(1e-3, 30), 'Blues'),
        ('Grüneisen parameter γ', o.gamma_Gruneisen, Normalize(0, 2.0), 'Blues'),
    ]
    fig, axs = plt.subplots(2, 3, figsize=(16, 9.6), sharex=True, sharey=True)
    for ax, (title, Z2, norm, cmap) in zip(axs.ravel(), panels):
        m = ax.pcolormesh(P, T, np.ma.masked_invalid(Z2).T, norm=norm, cmap=cmap, shading='auto',
                          rasterized=True)
        cb = fig.colorbar(m, ax=ax, pad=0.02, fraction=0.05)
        cb.ax.tick_params(labelsize=8)
        ok = np.isfinite(Tm)
        ax.fill_between(Pm[ok], 230, np.maximum(Tm[ok], 230), facecolor='none', hatch='////',
                        edgecolor='#8c8b87', lw=0)
        ax.plot(Pm[ok], Tm[ok], '-', color=INK, lw=1.0)
        ax.plot(sat.P, Tsat, '-', color='#eb6834', lw=1.8)
        ax.plot(22.064, 647.096, 'o', color='#eb6834', ms=5, mec=INK, mew=0.6)
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_xlim(P[0], P[-1]); ax.set_ylim(T[0], T[-1])
        ax.set_title(title, loc='left', color=INK, fontsize=10)
    for ax in axs[-1]:
        ax.set_xlabel('Pressure (MPa)')
    for ax in axs[:, 0]:
        ax.set_ylabel('Temperature (K)')
    ax0 = axs[0, 0]
    ax0.text(3e-6, 280, 'vapour', color=INK, fontsize=9)
    ax0.text(2e1, 300, 'liquid', color='white', fontsize=9)
    ax0.text(3e3, 250, 'solid\n(hatched)', color=INK, fontsize=8.5)
    fig.suptitle('water3 (psi surface stage5_23c) over its full range — stable branch; orange: saturation curve '
                 'and critical point; black: sublimation + melting curve of the stable solid', x=0.01, ha='left',
                 color=INK, fontsize=11)
    fig.tight_layout()
    f = os.path.join(outdir, 'water3_surface_maps.png')
    fig.savefig(f, dpi=140); plt.close(fig)
    return f


if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else '.'
    os.makedirs(out, exist_ok=True)
    print(make(out))
