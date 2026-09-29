"""SeaFreeze tutorial: the Helmholtz fluid 'water3', density input, and
liquid-vapour-ice coexistence.

    python3 water3_tutorial.py [OUTDIR]

Runs every example of the "water3" section of the Python README and saves
water3_tutorial_python.png in OUTDIR (default: this folder).
"""
import os
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import seafreeze as sf


def scatter(*cols):
    """1-D object array of tuples, the SeaFreeze scatter-input format."""
    pts = np.empty(len(cols[0]), dtype=object)
    for i, row in enumerate(zip(*cols)):
        pts[i] = tuple(float(v) for v in row)
    return pts


# ---------------------------------------------------------------------------
# 1. water3 at (P, T): same call and outputs as every other SeaFreeze phase
# ---------------------------------------------------------------------------
out = sf.getProp(scatter([0.101325, 100, 1000], [298.15, 300, 350]), 'water3')
print('1. water3 at 0.1/100/1000 MPa:')
print('   rho =', np.round(out.rho, 3), 'kg/m3;  Cp =', np.round(out.Cp, 1), 'J/kg/K;  vel =',
      np.round(out.vel, 1), 'm/s')

# grid input: rows are pressures, columns temperatures
grid = np.array([np.array([0.1, 50, 500]), np.array([280., 300, 320])], dtype=object)
g3 = sf.getProp(grid, 'water3', sf.seafreeze.defpath, 'rho', 'alpha')
g1 = sf.getProp(grid, 'water1', sf.seafreeze.defpath, 'rho', 'alpha')
print('   grid rho water3 - water1 (kg/m3):\n', np.round(g3.rho - g1.rho, 3))

# ---------------------------------------------------------------------------
# 2. One fluid, two branches: the stable phase is the lower Gibbs energy
# ---------------------------------------------------------------------------
T = np.linspace(280, 500, 221)
stable = sf.getProp(scatter(np.full(T.size, 0.101325), T), 'water3', sf.seafreeze.defpath, 'rho', 'G')
liquid = sf.getProp(scatter(np.full(T.size, 0.101325), T), 'water3', sf.seafreeze.defpath, 'rho', 'G',
                    branch='liquid')           # superheated liquid above 373.124 K
Tb = T[np.argmax(stable.rho < 100)]
print(f'2. at 0.101325 MPa the stable branch switches to vapour at {Tb:.0f} K '
      f'(saturation temperature 373.124 K)')

# ---------------------------------------------------------------------------
# 3. Density-temperature input, for water3 AND the Gibbs splines
# ---------------------------------------------------------------------------
rT = scatter([1000.0, 1100.0], [300.0, 300.0])
w3 = sf.getProp(rT, 'water3', sf.seafreeze.defpath, 'P', 'Cp', rhoT=True)       # direct
w1 = sf.getProp(rT, 'water1', sf.seafreeze.defpath, 'P', 'Cp', rhoT=True)       # via rho2P
print('3. P at 1000 and 1100 kg/m3, 300 K:  water3', np.round(w3.P, 3), 'MPa;  water1',
      np.round(w1.P, 3), 'MPa')
ice = sf.getProp(scatter([1330.0], [260.0]), 'VI', sf.seafreeze.defpath, 'P', 'Vp', 'Vs', rhoT=True)
brine = sf.getProp(scatter([1050.0], [300.0], [1.0]), 'NaClaq', sf.seafreeze.defpath, 'P', 'muw', rhoT=True)
print(f'   ice VI at 1330 kg/m3, 260 K: P = {ice.P[0]:.1f} MPa, Vp = {ice.Vp[0]:.0f} m/s;  '
      f'NaCl(aq) 1 mol/kg at 1050 kg/m3, 300 K: P = {brine.P[0]:.2f} MPa')
print('   rho2P (water3):', np.round(sf.rho2P([997.047, 1100], [298.15, 300], 'water3'), 4), 'MPa')

# ---------------------------------------------------------------------------
# 4. Phase equilibria with water3 as the liquid
# ---------------------------------------------------------------------------
pg = np.array([np.array([0.1, 300, 800]), np.array([260., 270, 276, 300])], dtype=object)
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    print('4. whichphase (0 = liquid, 1 = Ih, 6 = VI):\n', sf.whichphase(pg, 'water3'))
# down to the triple-point pressure (the default grid starts at 0.1 MPa)
line = sf.phase_lines('Ih', 'water3', P=np.geomspace(6.2e-4, 209, 300), T=np.arange(250, 273.3, 0.02))

# ---------------------------------------------------------------------------
# 5. Saturation and sublimation curves
# ---------------------------------------------------------------------------
sat = sf.saturation(np.linspace(230, 646.5, 200))       # < 273.16 K: supercooled liquid
Tsub = np.linspace(170, 273.16, 150)
sub = sf.sublimation(Tsub)                              # below 230 K: dilute-vapour extension (warns once)
print(f'5. p_sat(373.124 K) = {sf.saturation([373.124]).P[0]:.6f} MPa;  '
      f'p_sub(250 K) = {sf.sublimation([250.0]).P[0] * 1e6:.3f} Pa')

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
C3, C1, CV, INK2 = '#2a78d6', '#eb6834', '#1baf7a', '#52514e'
fig, axs = plt.subplots(2, 2, figsize=(12, 9))
ax = axs[0, 0]
ax.semilogy(T, stable.rho, '-', color=C3, lw=2.2, label="branch='stable' (default)")
ax.semilogy(T, liquid.rho, '--', color=C1, lw=1.6, label="branch='liquid' (metastable above T$_b$)")
ax.axvline(373.124, color=INK2, lw=0.8, ls='-.')
ax.set_xlabel('Temperature (K)'); ax.set_ylabel('Density (kg/m³)')
ax.set_title('Boiling at 0.101325 MPa: one EOS, two branches', loc='left')
ax.legend(fontsize=8.5)

ax = axs[0, 1]
rho = np.linspace(950, 1250, 61)
for Tq, ls in [(280.0, '-'), (350.0, '--')]:
    p3 = sf.getProp(np.array([rho, np.array([Tq])], dtype=object), 'water3', sf.seafreeze.defpath,
                    'P', rhoT=True).P[:, 0]
    p1 = sf.getProp(np.array([rho, np.array([Tq])], dtype=object), 'water1', sf.seafreeze.defpath,
                    'P', rhoT=True).P[:, 0]
    ax.plot(rho, p3, ls, color=C3, lw=2, label=f'water3, {Tq:.0f} K')
    ax.plot(rho, p1, ls, color=C1, lw=1.4, label=f'water1 (Gibbs, via rho2P), {Tq:.0f} K')
ax.set_xlabel('Density (kg/m³)'); ax.set_ylabel('Pressure (MPa)')
ax.set_title('(ρ,T) input: isotherms P(ρ)', loc='left')
ax.legend(fontsize=8)

ax = axs[1, 0]
ax.semilogy(sat.T, sat.P * 1e6, '-', color=C3, lw=2, label='liquid–vapour (water3)')
ax.semilogy(Tsub[Tsub >= 230], sub.P[Tsub >= 230] * 1e6, '-', color=CV, lw=2, label='ice Ih–vapour')
ax.semilogy(Tsub[Tsub < 230], sub.P[Tsub < 230] * 1e6, '--', color=CV, lw=1.6,
            label='ice Ih–vapour, dilute-vapour extension')
o = np.argsort(line.P)
ax.semilogy(line.T[o], line.P[o] * 1e6, '-', color=C1, lw=2, label='ice Ih–liquid (water3)')
ax.plot(273.16, 611.657, 'ks', ms=5)
ax.set_xlabel('Temperature (K)'); ax.set_ylabel('Pressure (Pa)')
ax.set_title('Triple point region: melting, boiling, sublimation', loc='left')
ax.set_ylim(1e-3, 3e8)
ax.legend(fontsize=8, loc='lower right')

ax = axs[1, 1]
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    sf.wpd(ax=ax, liquid='water3', phase_labels=True)
ax.set_title('')
ax.set_title("sf.wpd(liquid='water3')", loc='left')
ax.xaxis.label.set_size(10); ax.yaxis.label.set_size(10)
fig.tight_layout()

outdir = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
os.makedirs(outdir, exist_ok=True)
f = os.path.join(outdir, 'water3_tutorial_python.png')
fig.savefig(f, dpi=140)
print('saved', f)
