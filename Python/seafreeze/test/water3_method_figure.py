"""How SeaFreeze evaluates a Helmholtz fluid at (P,T): a four-panel explainer.

    python3 -m seafreeze.test.water3_method_figure OUTDIR

  A  one isotherm P(rho): the target P crosses the mechanically stable branches
     (dP/drho > 0) at a vapour-like and a liquid-like root
  B  G of both roots vs P: the lower one is the stable phase ('stable' branch);
     they cross at the saturation pressure
  C  near Tc the dome carries small (dP/drho)_T loops; a coarse bracket grid can
     land Newton on a loop root (the bug fixed by the fine dense-fluid grid)
  D  sublimation: Newton in ln P on G_vapour - G_ice = 0 at fixed T
"""
import os
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import seafreeze as sf
from seafreeze import coexistence as cx
from seafreeze.seafreeze import defpath, _load_spline
from lbftd import evalHelmholtz as eh

C3, CL, CV, CBAD = '#2a78d6', '#1d5aa6', '#eb6834', '#e34948'
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.grid': True,
                     'grid.color': GRID, 'grid.linewidth': 0.6, 'axes.spines.top': False,
                     'axes.spines.right': False, 'figure.facecolor': '#fcfcfb',
                     'axes.facecolor': '#fcfcfb', 'legend.frameon': False})


def _sc(P, T):
    out = np.empty(len(P), dtype=object)
    for i in range(len(P)):
        out[i] = (float(P[i]), float(T[i]))
    return out


def _isotherm(sp, rho, T):
    o = eh.evalHelmholtzScatter(sp, _sc(rho, np.full(rho.size, T)), 'P', 'Kt', 'G', rhoT=True)
    return o.P, o.Kt, o.G


def make(outdir):
    sp = _load_spline(defpath, 'water3')
    fig, axs = plt.subplots(2, 2, figsize=(13.5, 10))
    (axA, axB), (axC, axD) = axs

    # ---- A: isotherm at 400 K, target P = 0.1 MPa ------------------------------
    T, Pt = 400.0, 0.1
    rho = np.geomspace(1e-3, 1200, 4000)
    P, Kt, _ = _isotherm(sp, rho, T)
    stab = Kt > 0
    Pm = np.where(stab, P, np.nan); Pu = np.where(~stab, P, np.nan)
    axA.plot(rho, Pm, '-', color=C3, lw=2, label='dP/dρ > 0 (mechanically stable)')
    axA.plot(rho, Pu, '-', color=INK2, lw=1.2, alpha=0.7, label='dP/dρ < 0 (spinodal region)')
    axA.axhline(Pt, color=CV, lw=1.2, ls='--')
    rv = sf.getProp(_sc([Pt], [T]), 'water3', defpath, 'rho', branch='vapor').rho[0]
    rl = sf.getProp(_sc([Pt], [T]), 'water3', defpath, 'rho', branch='liquid').rho[0]
    axA.plot([rv], [Pt], 'o', color=CV, ms=9, mfc='white', mew=2, label=f'vapour root  ρ = {rv:.3f}')
    axA.plot([rl], [Pt], 's', color=CL, ms=9, mfc='white', mew=2, label=f'liquid root  ρ = {rl:.1f}')
    axA.set_xscale('log'); axA.set_yscale('symlog', linthresh=1e-2)
    axA.set_ylim(-500, 1e3)
    axA.set_xlabel('Density ρ (kg/m³)'); axA.set_ylabel('P(ρ, T) = ρ² ∂F/∂ρ  (MPa)')
    axA.set_title(f'A. Solve P(ρ,T) = {Pt} MPa at T = {T:.0f} K: bracket every upward crossing',
                  loc='left', color=INK, fontsize=10)
    axA.text(Pt and 2e-3, 0.13, 'target P', color=CV, fontsize=9)
    axA.legend(loc='lower right', fontsize=8.5)

    # ---- B: G of both roots vs P -----------------------------------------------
    Pg = np.geomspace(0.02, 1.0, 120)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        gl = sf.getProp(_sc(Pg, np.full(Pg.size, T)), 'water3', defpath, 'G', branch='liquid').G
        gv = sf.getProp(_sc(Pg, np.full(Pg.size, T)), 'water3', defpath, 'G', branch='vapor').G
    ps = sf.saturation([T]).P[0]
    axB.plot(Pg, gl / 1e3, '-', color=CL, lw=2, label='liquid root')
    axB.plot(Pg, gv / 1e3, '-', color=CV, lw=2, label='vapour root')
    axB.plot(Pg, np.minimum(gl, gv) / 1e3, ':', color=INK, lw=3, alpha=0.45, label="'stable' = min G")
    axB.axvline(ps, color=INK2, lw=1, ls='-.')
    axB.annotate(f'p$_σ$({T:.0f} K) = {ps * 1e3:.2f} kPa\n(IAPWS-95 aux: {cx.psat_iapws_aux(T) * 1e3:.2f})',
                 (ps, np.interp(ps, Pg, gl) / 1e3), xytext=(12, -40), textcoords='offset points',
                 color=INK, fontsize=9, arrowprops=dict(arrowstyle='-', color=INK2))
    axB.set_xscale('log')
    axB.set_xlabel('Pressure (MPa)'); axB.set_ylabel('G (kJ/kg)')
    axB.set_title('B. Which root? The lower Gibbs energy is the stable phase', loc='left', color=INK, fontsize=10)
    axB.legend(loc='upper left', fontsize=8.5)

    # ---- C: near-critical loops and the bracket grid ---------------------------
    T2 = 634.5
    rho2 = np.linspace(80, 620, 5000)
    P2, Kt2, _ = _isotherm(sp, rho2, T2)
    s2 = sf.saturation([T2])
    axC.plot(rho2, P2, '-', color=C3, lw=1.8, label=f'water3 isotherm {T2} K')
    axC.axhline(s2.P[0], color=INK2, lw=1, ls='-.', label=f'p$_σ$ = {s2.P[0]:.3f} MPa (Maxwell)')
    xs = np.geomspace(1e-12, 16000, 600)
    xc = xs[(xs > 80) & (xs < 620)]
    Pc, _, _ = _isotherm(sp, xc, T2)
    axC.plot(xc, Pc, 'x', color=CBAD, ms=7, mew=1.6, label='coarse bracket grid (6 % steps): misses loops')
    axC.plot(s2.rho_A, s2.P, 's', color=CL, ms=9, mfc='white', mew=2, label=f'true liquid root ρ = {s2.rho_A[0]:.1f}')
    axC.plot(s2.rho_B, s2.P, 'o', color=CV, ms=9, mfc='white', mew=2, label=f'true vapour root ρ = {s2.rho_B[0]:.1f}')
    axC.plot([486.584], [18.65620], 'v', color=CBAD, ms=10, label='old result: loop root ρ = 486.6 (p$_σ$ −1.7 %)')
    axC.set_ylim(s2.P[0] - 0.6, s2.P[0] + 0.6)
    axC.set_xlabel('Density ρ (kg/m³)'); axC.set_ylabel('P (MPa)')
    axC.set_title('C. Near T$_c$: small (∂P/∂ρ)$_T$ loops in the dome → fine 2 kg/m³ bracket grid',
                  loc='left', color=INK, fontsize=10)
    axC.legend(loc='upper left', fontsize=8, framealpha=0.92, frameon=True, facecolor='#fcfcfb', edgecolor='none')
    axC.text(0.99, 0.02, 'dome interior is data-free in the fit:\nloops are expected there, not physical', transform=axC.transAxes,
             ha='right', va='bottom', fontsize=8.5, color=INK2)

    # ---- D: sublimation by Newton in ln P ---------------------------------------
    T3 = 250.0
    Pd = np.geomspace(20e-6, 300e-6, 100)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        gi = sf.getProp(_sc(Pd, np.full(Pd.size, T3)), 'Ih', defpath, 'G').G
        gvv = sf.getProp(_sc(Pd, np.full(Pd.size, T3)), 'water3', defpath, 'G', branch='vapor').G
    axD.plot(Pd * 1e6, (gvv - gi) / 1e3, '-', color=C3, lw=2, label='G$_{vapour}$(water3) − G$_{ice Ih}$')
    axD.axhline(0, color=INK2, lw=0.8)
    # Newton iterates from a deliberately poor start (4x R14-08)
    lnP = np.log(4 * cx.psub_iapws(T3)); its = []
    for _ in range(6):
        P = np.exp(lnP)
        a = sf.getProp(_sc([P], [T3]), 'Ih', defpath, 'G', 'rho')
        b = sf.getProp(_sc([P], [T3]), 'water3', defpath, 'G', 'rho', branch='vapor')
        its.append((P, (b.G[0] - a.G[0])))
        lnP -= (a.G[0] - b.G[0]) / (P * 1e6 * (1 / a.rho[0] - 1 / b.rho[0]))
    for k, (P, dg) in enumerate(its):
        axD.plot(P * 1e6, dg / 1e3, 'o', color=CV, ms=7 if k else 9, mfc='white', mew=1.6)
        if k == 0:
            axD.annotate('start (4 × R14-08)', (P * 1e6, dg / 1e3), xytext=(-110, -4), textcoords='offset points',
                         color=INK2, fontsize=8.5)
    ps3 = sf.sublimation([T3]).P[0]
    axD.axvline(ps3 * 1e6, color=INK2, lw=1, ls='-.')
    axD.text(0.03, 0.97, f'p$_{{sub}}$({T3:.0f} K) = {ps3 * 1e6:.3f} Pa  (R14-08: {cx.psub_iapws(T3) * 1e6:.3f} Pa)\n'
             f'converged in {sum(abs(d) > 1e-6 for _, d in its)} Newton steps',
             transform=axD.transAxes, va='top', color=INK, fontsize=9)
    axD.set_xscale('log')
    axD.set_xlabel('Pressure (Pa)'); axD.set_ylabel('ΔG (kJ/kg)')
    axD.set_title('D. Sublimation: Newton in ln P on G$_{vap}$ = G$_{Ih}$ (slope = P(1/ρ$_{ice}$ − 1/ρ$_{vap}$))',
                  loc='left', color=INK, fontsize=10)
    axD.legend(loc='lower right', fontsize=8.5)

    fig.tight_layout()
    f = os.path.join(outdir, 'water3_method.png')
    fig.savefig(f, dpi=150); plt.close(fig)
    return f


if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else '.'
    os.makedirs(out, exist_ok=True)
    print(make(out))
