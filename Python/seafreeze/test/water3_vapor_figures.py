"""Vapour-side validation figures for the Helmholtz fluid 'water_Brown2026'.

    python3 -m seafreeze.test.water3_vapor_figures OUTDIR [NISTDIR] [IAPWS95_DATA_DIR]

Writes to OUTDIR:
  water3_sublimation.png     ice Ih sublimation pressure: water_Brown2026 (+ dilute-vapour
                             extension below 230 K) vs IAPWS R14-08 and the NIST
                             measurements of Bielska et al. (2013); supercooled-
                             liquid vapour pressure for context
  water3_vapor_density.png   saturated-vapour density and vapour isotherms vs the
                             NIST Chemistry WebBook (IAPWS-95) and the saturated-
                             vapour density measurements used to build IAPWS-95

NISTDIR caches NIST WebBook tables (fetched on first use).  IAPWS95_DATA_DIR is
lbf-thermo's materials/H2O_liquid/data/IAPWS95_source_data (normalized.csv);
the experimental overlay is skipped without it.
"""
import os
import sys
import urllib.request
import warnings

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import seafreeze as sf
from seafreeze import coexistence as cx
from seafreeze.seafreeze import defpath

C3, CREF, CEXP, CDATA = '#2a78d6', '#eb6834', '#1baf7a', '#4a3aa7'
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.grid': True,
                     'grid.color': GRID, 'grid.linewidth': 0.6, 'axes.spines.top': False,
                     'axes.spines.right': False, 'figure.facecolor': '#fcfcfb',
                     'axes.facecolor': '#fcfcfb', 'legend.frameon': False})

# Bielska, Havey, Scace, Lisak, Harvey & Hodges (2013), Geophys. Res. Lett. 40,
# 6303-6307, doi:10.1002/2013GL058474 (NIST), Table 2: ice Ih vapour pressure,
# referenced to p_t = 611.657 Pa.  T/K, p/Pa, u(p)/Pa (k = 1).
BIELSKA2013 = np.array([
    [253.380, 1.055e2, 4.1e-1], [238.375, 2.289e1, 8.9e-2], [223.521, 4.134e0, 1.7e-2],
    [204.937, 3.395e-1, 1.4e-3], [204.918, 3.388e-1, 1.5e-3], [204.918, 3.398e-1, 1.5e-3],
    [201.230, 1.959e-1, 8.4e-4], [196.757, 9.808e-2, 4.2e-4], [192.476, 4.892e-2, 2.1e-4],
    [187.065, 1.952e-2, 8.5e-5], [183.205, 9.766e-3, 4.5e-5], [179.478, 4.912e-3, 2.7e-5],
    [174.777, 1.948e-3, 1.4e-5]])

WEBBOOK = ('https://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID=C7732185&Digits=8'
           '&RefState=DEF&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs'
           '&VisUnit=uPa*s&STUnit=N%2Fm')
ISOTHERMS = {300: (0.0001, 0.0035, 0.0001), 400: (0.001, 0.24, 0.005), 500: (0.01, 2.6, 0.05),
             600: (0.1, 12.3, 0.25), 700: (0.1, 60, 1), 1000: (0.1, 100, 2)}


def _webbook(nistdir, name, query):
    f = os.path.join(nistdir, name)
    if not os.path.exists(f):
        os.makedirs(nistdir, exist_ok=True)
        with urllib.request.urlopen(WEBBOOK + query, timeout=60) as r, open(f, 'wb') as out:
            out.write(r.read())
    with open(f) as fh:
        head = fh.readline().rstrip('\n').split('\t')
        rows = [ln.rstrip('\n').split('\t') for ln in fh if ln.strip()]
    return head, rows


def _col(head, rows, name, as_float=True):
    j = head.index(name)
    return np.array([float(r[j]) for r in rows]) if as_float else np.array([r[j] for r in rows])


def _scatter(P, T):
    out = np.empty(len(P), dtype=object)
    for i in range(len(P)):
        out[i] = (float(P[i]), float(T[i]))
    return out


def sublimation_figure(outdir):
    Tfull = np.linspace(230.0, 273.16, 88)
    Text = np.linspace(170.0, 230.0, 121)
    full = sf.sublimation(Tfull)
    ext = sf.sublimation(Text, dilute_extension=True)
    Tl = np.linspace(230.0, 300.0, 71)
    sat = sf.saturation(Tl)                                   # < 273.16 K: supercooled liquid
    Tb, pb, ub = BIELSKA2013.T
    pref_b = cx.psub_iapws(Tb) * 1e6

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))
    ax.semilogy(Text, cx.psub_iapws(Text) * 1e6, ':', color=CREF, lw=1.6)
    ax.semilogy(Tfull, cx.psub_iapws(Tfull) * 1e6, ':', color=CREF, lw=1.6, label='IAPWS R14-08 (Wagner et al. 2011)')
    ax.semilogy(Tfull, full.P * 1e6, '-', color=C3, lw=2.2, label='water_Brown2026 vapour + ice Ih (SeaFreeze)')
    ax.semilogy(Text, ext.P * 1e6, '--', color=C3, lw=1.6, label='water_Brown2026 dilute-vapour extension (< 230 K)')
    ax.semilogy(Tl, sat.P * 1e6, '-', color=INK2, lw=1.2, label='water_Brown2026 liquid–vapour (supercooled < 273.16 K)')
    ax.errorbar(Tb, pb, yerr=ub, fmt='o', ms=5, color=CDATA, mfc='white', mew=1.4,
                label='NIST measurements (Bielska et al. 2013)')
    ax.plot(273.16, 611.657, 's', color=INK, ms=6)
    ax.annotate('triple point', (273.16, 611.657), xytext=(-70, 8), textcoords='offset points', color=INK2)
    ax.axvline(230, color=INK2, lw=0.8, ls='-.')
    ax.text(231, 3e2, 'water_Brown2026\nT$_{min}$ = 230 K', color=INK2, fontsize=9, va='top')
    ax.set_xlabel('Temperature (K)'); ax.set_ylabel('Vapour pressure (Pa)')
    ax.set_title('Sublimation pressure of ice Ih', loc='left', color=INK)
    ax.legend(loc='lower right', fontsize=8.5)
    ax.set_xlim(170, 302)

    ax2.axhline(0, color=CREF, lw=1.2, ls=':')
    ax2.plot(Tfull, 100 * (full.P / cx.psub_iapws(Tfull) - 1), '-', color=C3, lw=2.2, label='water_Brown2026 (full surface)')
    ax2.plot(Text, 100 * (ext.P / cx.psub_iapws(Text) - 1), '--', color=C3, lw=1.6, label='water_Brown2026 dilute-vapour extension')
    ax2.errorbar(Tb, 100 * (pb / pref_b - 1), yerr=100 * ub / pref_b, fmt='o', ms=5, color=CDATA,
                 mfc='white', mew=1.4, capsize=2, label='NIST Bielska et al. 2013 (±1σ)')
    ax2.axvline(230, color=INK2, lw=0.8, ls='-.')
    ax2.set_xlabel('Temperature (K)'); ax2.set_ylabel('100 · (p / p$_{R14-08}$ − 1)  (%)')
    ax2.set_title('Deviation from IAPWS R14-08', loc='left', color=INK)
    ax2.legend(loc='upper left', fontsize=8.5)
    ax2.set_ylim(-1.0, 1.0)
    fig.tight_layout()
    f = os.path.join(outdir, 'water3_sublimation.png')
    fig.savefig(f, dpi=150); plt.close(fig)
    dev = 100 * (full.P / cx.psub_iapws(Tfull) - 1)
    print(f'sublimation 230-273.16 K: water_Brown2026 vs R14-08 {dev.min():+.4f} .. {dev.max():+.4f} %')
    return f


def vapor_density_figure(outdir, nistdir, datadir):
    head, rows = _webbook(nistdir, 'satT.txt', '&Type=SatT&TLow=275&THigh=646&TInc=5')
    Tn = _col(head, rows, 'Temperature (K)')
    Pn = _col(head, rows, 'Pressure (MPa)')
    rvn = _col(head, rows, 'Density (v, kg/m3)')
    o = np.argsort(Tn); Tn, Pn, rvn = Tn[o], Pn[o], rvn[o]
    sat = sf.saturation(Tn)
    Tf = np.linspace(273.16, 646.5, 300)
    satf = sf.saturation(Tf)

    fig, axs = plt.subplots(2, 2, figsize=(13, 9.5))
    ax, ax2, ax3, ax4 = axs.ravel()
    ax.semilogy(Tf, satf.rho_B, '-', color=C3, lw=2.2, label='water_Brown2026 (SeaFreeze)')
    ax.semilogy(Tn[::12], rvn[::12], 'o', ms=4, color=CREF, mfc='white', mew=1.2,
                label='NIST WebBook (IAPWS-95)')
    ax2.plot(Tn, 100 * (sat.rho_B / rvn - 1), 'o', ms=3, color=CREF, mfc='white', mew=1.0,
             label='ρ$_v$ vs NIST WebBook (IAPWS-95)')
    ax2.plot(Tn, 100 * (sat.P / Pn - 1), '.', ms=3, color=INK2, label='p$_\\sigma$ vs NIST WebBook')
    if datadir and os.path.exists(os.path.join(datadir, 'normalized.csv')):
        import pandas as pd
        d = pd.read_csv(os.path.join(datadir, 'normalized.csv'))
        dv = d[(d.property == 'rhosat_v') & (d.T_K < 646)]
        se = sf.saturation(dv.T_K.values)
        ax.semilogy(dv.T_K, dv.value, '^', ms=4, color=CEXP, mfc='none', mew=1.0,
                    label='Osborne et al. 1937/39 (IAPWS-95 source data)')
        ax2.plot(dv.T_K, 100 * (se.rho_B / dv.value.values - 1), '^', ms=4, color=CEXP, mfc='none',
                 mew=1.0, label='vs Osborne et al. 1937/39 (measured)')
    ax.set_xlabel('Temperature (K)'); ax.set_ylabel('Saturated vapour density (kg/m³)')
    ax.set_title('Saturated vapour density', loc='left', color=INK)
    ax.legend(loc='lower right', fontsize=8.5)
    ax2.axhline(0, color=INK2, lw=0.8)
    ax2.set_xlabel('Temperature (K)'); ax2.set_ylabel('100 · (water_Brown2026 / reference − 1)  (%)')
    ax2.set_title('Saturated vapour: deviations', loc='left', color=INK)
    ax2.set_ylim(-0.3, 0.6)
    ax2.legend(loc='lower left', fontsize=8.5)

    shades = ['#b7d3f4', '#8ab8ee', '#5d9de6', '#2a78d6', '#1d5aa6', '#123c70']   # one hue, light -> dark
    for (T, (lo, hi, inc)), cshade in zip(ISOTHERMS.items(), shades):
        head, rows = _webbook(nistdir, f'iso_{T}.txt',
                              f'&Type=IsoTherm&T={T}&PLow={lo}&PHigh={hi}&PInc={inc}')
        ph = _col(head, rows, 'Phase', as_float=False)
        keep = np.isin(ph, ['vapor', 'supercritical'])
        P = _col(head, rows, 'Pressure (MPa)')[keep]
        rn = _col(head, rows, 'Density (kg/m3)')[keep]
        o = np.argsort(P); P, rn = P[o], rn[o]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            r3 = sf.getProp(_scatter(P, np.full(P.size, T)), 'water_Brown2026', defpath, 'rho', branch='vapor').rho
        ax3.loglog(P, r3, '-', color=cshade, lw=2)
        ax3.loglog(P[::4], rn[::4], 'o', ms=3.5, color=CREF, mfc='white', mew=1.0)
        ax3.annotate(f'{T} K', (P[-1], r3[-1]), xytext=(4, 0), textcoords='offset points',
                     color=INK2, fontsize=9, va='center')
        ax4.semilogx(P, 100 * (r3 / rn - 1), '-', color=cshade, lw=2, label=f'{T} K')
    ax3.plot([], [], '-', color=C3, lw=2, label='water_Brown2026 (SeaFreeze)')
    ax3.plot([], [], 'o', ms=3.5, color=CREF, mfc='white', mew=1.0, label='NIST WebBook (IAPWS-95)')
    ax3.set_xlabel('Pressure (MPa)'); ax3.set_ylabel('Density (kg/m³)')
    ax3.set_title('Vapour / supercritical isotherms (to saturation)', loc='left', color=INK)
    ax3.legend(loc='upper left', fontsize=8.5)
    ax4.axhline(0, color=INK2, lw=0.8)
    ax4.set_xlabel('Pressure (MPa)'); ax4.set_ylabel('100 · (ρ$_{water_Brown2026}$ / ρ$_{NIST}$ − 1)  (%)')
    ax4.set_title('Isotherm density deviations from NIST WebBook', loc='left', color=INK)
    ax4.legend(loc='upper left', fontsize=8.5, ncol=2)
    fig.tight_layout()
    f = os.path.join(outdir, 'water3_vapor_density.png')
    fig.savefig(f, dpi=150); plt.close(fig)
    d = 100 * (sat.rho_B / rvn - 1)
    print(f'saturated vapour density vs NIST: {np.nanmin(d):+.4f} .. {np.nanmax(d):+.4f} % (275-646 K)')
    return f


if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else '.'
    nist = sys.argv[2] if len(sys.argv) > 2 else os.path.join(out, 'nist_webbook')
    data = sys.argv[3] if len(sys.argv) > 3 else None
    os.makedirs(out, exist_ok=True)
    print(sublimation_figure(out))
    print(vapor_density_figure(out, nist, data))
