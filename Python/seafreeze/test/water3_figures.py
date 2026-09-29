"""Diagnostic figures for the Helmholtz liquid 'water3' (psi surface).

    python3 -m seafreeze.test.water3_figures [outdir]

Writes, to outdir (default: current directory):
  water3_phase_diagram.png   ice–liquid melting curves, water3 vs water1, plus
                             the melting-temperature difference along each curve
  water3_properties.png      rho, Cp, sound speed, alpha along isobars for
                             water3, water1 and water_IAPWS95 (Gibbs spline)
  water3_deviation_map.png   relative difference water3 - water1 in rho and Cp
                             over the water1 domain
"""
import os
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import seafreeze as sf

# Categorical slots (fixed order): water3, water1, IAPWS-95 — validated palette
C = {'water3': '#2a78d6', 'water1': '#eb6834', 'water_IAPWS95': '#1baf7a'}
LS = {'water3': '-', 'water1': '--', 'water_IAPWS95': ':'}
LABEL = {'water3': 'water3 (Helmholtz psi, 2026)', 'water1': 'water1 (Bollengier 2019)',
         'water_IAPWS95': 'IAPWS-95 (Gibbs spline)'}
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.grid': True,
                     'grid.color': GRID, 'grid.linewidth': 0.6, 'axes.spines.top': False,
                     'axes.spines.right': False, 'figure.facecolor': '#fcfcfb',
                     'axes.facecolor': '#fcfcfb', 'legend.frameon': False})


def grid(P, T):
    return np.array([np.asarray(P, float), np.asarray(T, float)], dtype=object)


def phase_diagram(outdir):
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 5.2), gridspec_kw={'width_ratios': [1.35, 1]})
    # ice–ice boundaries (liquid-independent)
    for a, b in [('Ih', 'II'), ('Ih', 'III'), ('II', 'III'), ('II', 'V'), ('II', 'VI'), ('III', 'V'), ('V', 'VI')]:
        r = sf.phase_lines(a, b, segment='stable')
        o = np.argsort(r.T)
        ax.plot(r.P[o], r.T[o], '.', color=INK2, ms=1.2)
    for ice in ['Ih', 'III', 'V', 'VI']:
        res = {}
        for liq in ['water1', 'water3']:
            r = sf.phase_lines(ice, liq, segment='stable')
            o = np.argsort(r.P)
            res[liq] = (r.P[o], r.T[o])
            ax.plot(r.P[o], r.T[o], LS[liq], color=C[liq], lw=2 if liq == 'water3' else 1.6,
                    label=LABEL[liq] if ice == 'Ih' else None)
        P3, T3 = res['water3']; P1, T1 = res['water1']
        Pc = np.linspace(max(P3.min(), P1.min()), min(P3.max(), P1.max()), 200)
        dT = np.interp(Pc, P3, T3) - np.interp(Pc, P1, T1)
        ax2.plot(Pc, 1e3 * dT, '-', color=C['water3'], lw=2)
        ax2.annotate(ice, (Pc[len(Pc) // 2], 1e3 * dT[len(Pc) // 2]), textcoords='offset points',
                     xytext=(0, 6), ha='center', color=INK2, fontsize=9)
    for txt, (p, t) in {'Ih': (80, 240), 'II': (330, 217), 'III': (280, 248), 'V': (490, 247),
                        'VI': (1200, 255), 'Liquid': (300, 330)}.items():
        ax.text(p, t, txt, color=INK, fontsize=10, fontweight='bold', ha='center')
    ax.set_xlim(0, 2300); ax.set_ylim(180, 380)
    ax.set_xlabel('Pressure (MPa)'); ax.set_ylabel('Temperature (K)')
    ax.set_title('Ice–liquid equilibrium: water3 vs water1', loc='left', color=INK)
    ax.legend(loc='lower right')
    ax2.axhline(0, color=INK2, lw=0.8)
    ax2.set_xlabel('Pressure (MPa)'); ax2.set_ylabel('T_melt(water3) − T_melt(water1)  (mK)')
    ax2.set_title('Melting-temperature difference (stable segments)', loc='left', color=INK)
    fig.tight_layout()
    f = os.path.join(outdir, 'water3_phase_diagram.png')
    fig.savefig(f, dpi=150); plt.close(fig)
    return f


def properties(outdir):
    isobars = [0.1, 400.0, 1000.0]
    props = [('rho', 'Density (kg/m³)'), ('Cp', 'Cp (J/kg/K)'), ('vel', 'Sound speed (m/s)'),
             ('alpha', 'α (1/K)')]
    T = np.linspace(240, 500, 131)
    fig, axs = plt.subplots(len(props), len(isobars), figsize=(12, 11), sharex=True)
    for j, P in enumerate(isobars):
        for mat in ['water_IAPWS95', 'water1', 'water3']:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                o = sf.getProp(grid([P], T), mat)
            for i, (k, lab) in enumerate(props):
                y = np.asarray(getattr(o, k)).ravel()
                axs[i, j].plot(T, y, LS[mat], color=C[mat], lw=2 if mat == 'water3' else 1.5,
                               label=LABEL[mat])
        axs[0, j].set_title(f'P = {P:g} MPa', loc='left', color=INK)
        for i, (k, lab) in enumerate(props):
            if j == 0:
                axs[i, j].set_ylabel(lab)
        axs[-1, j].set_xlabel('Temperature (K)')
    axs[0, 0].legend(loc='lower left', fontsize=8.5)
    # Cp at 0.1 MPa below ~260 K: water1 and IAPWS-95 are extrapolating
    fig.suptitle('Liquid-water properties along isobars', x=0.01, ha='left', color=INK)
    fig.tight_layout()
    f = os.path.join(outdir, 'water3_properties.png')
    fig.savefig(f, dpi=150); plt.close(fig)
    return f


def deviation_map(outdir):
    P = np.linspace(0.1, 2300, 116)
    T = np.linspace(240, 500, 105)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        w3 = sf.getProp(grid(P, T), 'water3', sf.seafreeze.defpath, 'rho', 'Cp')
        w1 = sf.getProp(grid(P, T), 'water1', sf.seafreeze.defpath, 'rho', 'Cp')
        ph = sf.whichphase(grid(P, T))
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    for ax, k, lim, lab in [(axs[0], 'rho', 0.2, 'ρ'), (axs[1], 'Cp', 3.0, 'Cp')]:
        d = 100 * (getattr(w3, k) / getattr(w1, k) - 1)
        d = np.where(ph == 0, d, np.nan)            # stable liquid only
        m = ax.pcolormesh(P, T, d.T, cmap='RdBu_r', vmin=-lim, vmax=lim, shading='auto')
        cb = fig.colorbar(m, ax=ax); cb.set_label(f'100·({lab}_water3/{lab}_water1 − 1)  (%)')
        ax.set_title(f'{lab}: water3 vs water1 (stable liquid)', loc='left', color=INK)
        ax.set_xlabel('Pressure (MPa)'); ax.grid(False)
        rms = np.sqrt(np.nanmean(d ** 2)); mx = np.nanmax(np.abs(d))
        ax.text(0.02, 0.97, f'rms {rms:.3g} %   max |Δ| {mx:.3g} %', transform=ax.transAxes,
                va='top', color=INK, fontsize=9, bbox=dict(fc='#fcfcfb', ec='none', alpha=0.8))
    axs[0].set_ylabel('Temperature (K)')
    fig.tight_layout()
    f = os.path.join(outdir, 'water3_deviation_map.png')
    fig.savefig(f, dpi=150); plt.close(fig)
    return f


if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else '.'
    os.makedirs(out, exist_ok=True)
    for fn in (phase_diagram, properties, deviation_map):
        print(fn(out))
