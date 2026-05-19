"""
Compare Python phase_lines vs Matlab reference.

Run from Python/ directory:
    python seafreeze/test/plot_phaselines_comparison.py
"""

import os
import sys
import warnings
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from seafreeze.phaselines import phase_lines, phase_range

REF_PATH = os.path.normpath(os.path.join(
    os.path.dirname(__file__), '..', '..', '..', 'Matlab', 'test', 'reference_phaselines.mat'))

ref = loadmat(REF_PATH, squeeze_me=True, struct_as_record=False)

# All case names in the reference
PURE_PAIRS = [
    'Ih_water1', 'Ih_II', 'Ih_III',
    'II_III', 'II_V', 'II_VI', 'II_water1',
    'III_V', 'III_water1',
    'V_water1', 'V_VI',
    'VI_water1',
]
NACL_PAIRS = [
    'Ih_NaClaq_m0p5', 'Ih_NaClaq_m2',
    'III_NaClaq_m0p5', 'III_NaClaq_m2',
    'V_NaClaq_m0p5', 'V_NaClaq_m2',
]

ALL_CASES = PURE_PAIRS + NACL_PAIRS

def compute_python_curve(case_name):
    case = ref[case_name]
    matA = str(case.matA)
    matB = str(case.matB)
    m_val = float(case.m) if not np.isnan(case.m) else None
    dT = float(case.dT)

    rA = phase_range(matA)
    rB = phase_range(matB)
    Tlo = max(rA.T[0], rB.T[0], 1.0)
    Thi = min(rA.T[1], rB.T[1])
    T_grid = np.arange(Tlo, Thi + dT / 2.0, dT)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        if m_val is not None:
            r = phase_lines(matA, matB, m=m_val, T=T_grid)
        else:
            r = phase_lines(matA, matB, T=T_grid)

    return (np.asarray(r.T, dtype=float).ravel(),
            np.asarray(r.P, dtype=float).ravel())


# ---- Figure 1: Pure-water pairs ----
fig1, axes1 = plt.subplots(3, 4, figsize=(16, 12))
fig1.suptitle('Pure-water phase lines: Matlab (blue) vs Python (red dashed)', fontsize=13)

for ax, name in zip(axes1.ravel(), PURE_PAIRS):
    case = ref[name]
    T_ref = np.asarray(case.T, dtype=float).ravel()
    P_ref = np.asarray(case.P, dtype=float).ravel()
    T_py, P_py = compute_python_curve(name)
    ax.plot(T_ref, P_ref, 'b-', lw=1.5, label='Matlab')
    ax.plot(T_py,  P_py,  'r--', lw=1.0, label='Python')
    ax.set_title(name.replace('_', ' '), fontsize=8)
    ax.set_xlabel('T (K)', fontsize=7)
    ax.set_ylabel('P (MPa)', fontsize=7)
    ax.tick_params(labelsize=6)

# hide spare axes
for ax in axes1.ravel()[len(PURE_PAIRS):]:
    ax.set_visible(False)

axes1.ravel()[0].legend(fontsize=7)
fig1.tight_layout()
out1 = os.path.join(os.path.dirname(__file__), 'phaselines_pure_comparison.png')
fig1.savefig(out1, dpi=150)
print(f"Saved {out1}")

# ---- Figure 2: NaClaq pairs ----
fig2, axes2 = plt.subplots(2, 3, figsize=(14, 9))
fig2.suptitle('NaClaq phase lines: Matlab (blue) vs Python (red dashed)', fontsize=13)

for ax, name in zip(axes2.ravel(), NACL_PAIRS):
    case = ref[name]
    T_ref = np.asarray(case.T, dtype=float).ravel()
    P_ref = np.asarray(case.P, dtype=float).ravel()
    T_py, P_py = compute_python_curve(name)
    ax.plot(T_ref, P_ref, 'b-', lw=1.5, label='Matlab')
    ax.plot(T_py,  P_py,  'r--', lw=1.0, label='Python')
    m_str = name.split('_m')[-1].replace('p', '.')
    ax.set_title(f"{name.split('_NaClaq')[0]} / NaClaq  m={m_str} mol/kg", fontsize=8)
    ax.set_xlabel('T (K)', fontsize=7)
    ax.set_ylabel('P (MPa)', fontsize=7)
    ax.tick_params(labelsize=6)

axes2.ravel()[0].legend(fontsize=7)
fig2.tight_layout()
out2 = os.path.join(os.path.dirname(__file__), 'phaselines_nacl_comparison.png')
fig2.savefig(out2, dpi=150)
print(f"Saved {out2}")
