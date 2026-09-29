"""
Generate a .mat reference file from the Python SeaFreeze implementation,
for consumption by Matlab/test/test_SeaFreeze_vs_python.m.

Baptiste Journaux - 2026

Run once from this directory (Matlab/test/):

    python gen_python_reference.py

Requires: numpy, scipy, hdf5storage, and the project's Python/ directory
resolvable (this script prepends it to sys.path automatically).

Notes on expected agreement with the MATLAB implementation:
  * For pure phases (Ih, II, III, V, VI, VII_X_French, water1, water_IAPWS95)
    the math is identical; any small numerical differences come from the
    spline evaluators (scipy splev vs. MATLAB fnval) and from any
    Bform->LBF representation mismatch between the two .mat files.
  * For NaClaq, MATLAB now prepends m=cutoff (2e-4) to match Python's
    convention, so apparent quantities (Cpa, Va, Vex, phi, aw) and the
    base/partial-molar properties (G, rho, Cp, mus, muw, Vm, Cpm) should
    all match. The MATLAB port no longer computes gamma/G_ss/Gex — those
    are still produced on the Python side but are not part of the
    cross-version comparison.
"""

import os
import sys
import numpy as np
from scipy.io import savemat

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'Python'))

import warnings  # noqa: E402
from seafreeze import seafreeze as sf  # noqa: E402


def direct_getProp(PTm, phase, tdv_names, rhoT=False):
    """Evaluate through the public sf.getProp, exactly as users call it
    (stitched NaClaq, Helmholtz water3, (rho,T) input).  The earlier bypass of
    getProp worked around a numpy-2 shear bug fixed in 1.1.2 and relied on a
    PhaseDesc field that no longer exists.
    """
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return sf.getProp(PTm, phase, sf.defpath, *tdv_names, rhoT=rhoT)


def grid(P, T, m=None):
    arr = [np.asarray(P, dtype=float), np.asarray(T, dtype=float)]
    if m is not None:
        arr.append(np.asarray(m, dtype=float))
    return np.array(arr, dtype=object)


def scatter(P, T, m=None):
    n = len(P)
    pts = np.empty(n, dtype=object)
    for i in range(n):
        if m is None:
            pts[i] = (float(P[i]), float(T[i]))
        else:
            pts[i] = (float(P[i]), float(T[i]), float(m[i]))
    return pts


BASE = ['G', 'S', 'U', 'H', 'A', 'rho', 'Cp', 'Cv', 'Kt', 'Kp', 'Ks', 'alpha', 'vel']
NACL = BASE + ['mus', 'muw', 'Vm', 'Cpm', 'Cpa', 'phi', 'aw', 'Va', 'Vex']


def collect(out, names):
    d = {}
    for n in names:
        if hasattr(out, n):
            arr = np.asarray(getattr(out, n))
            # When `m` includes the cutoff value, lbftd returns dtype=object
            # arrays (each cell a 1x1 ndarray). Flatten to plain float64 so
            # scipy.io.savemat writes a numeric matrix, not a cell array.
            if arr.dtype == object:
                arr = np.array(arr.tolist(), dtype=float).reshape(arr.shape)
            d[n] = arr
    return d


def add_case(cases, label, phase, PTm, names, extras, rhoT=False):
    out = direct_getProp(PTm, phase, names + (['P'] if rhoT else []), rhoT=rhoT)
    d = collect(out, names)
    d.update(extras)
    if rhoT:
        d['P_calc'] = np.asarray(out.P, dtype=float)   # computed pressure (MPa)
    d['phase'] = phase  # store phase explicitly so the MATLAB test
                        # doesn't have to parse the case name
    d['rhoT'] = int(rhoT)   # 1: the P field holds densities (kg/m3), P is an output
    cases[label] = d
    print(f"  + {label:24s} ({phase}, {len(names)} TDVs)")


def main():
    cases = {}
    print("Generating Python reference values...")

    # ---- Ice Ih ----
    P = np.arange(0.1, 200.1, 50.0)
    T = np.arange(240.0, 274.0, 10.0)
    add_case(cases, 'Ih_grid', 'Ih', grid(P, T), BASE, {'P': P, 'T': T})
    Ps = np.array([50.0, 100.0, 150.0])
    Ts = np.array([250.0, 255.0, 265.0])
    add_case(cases, 'Ih_scatter', 'Ih', scatter(Ps, Ts), BASE, {'P': Ps, 'T': Ts})

    # ---- Ice II, III, V, VI ----
    P = np.array([250.0, 300.0, 350.0])
    T = np.array([240.0, 245.0, 250.0])
    add_case(cases, 'II_grid', 'II', grid(P, T), BASE, {'P': P, 'T': T})

    P = np.array([220.0, 260.0, 340.0])
    T = np.array([240.0, 245.0, 250.0])
    add_case(cases, 'III_grid', 'III', grid(P, T), BASE, {'P': P, 'T': T})

    P = np.array([400.0, 500.0, 600.0])
    T = np.array([240.0, 250.0, 260.0])
    add_case(cases, 'V_grid', 'V', grid(P, T), BASE, {'P': P, 'T': T})

    P = np.array([700.0, 900.0, 1100.0])
    T = np.array([250.0, 270.0, 290.0])
    add_case(cases, 'VI_grid', 'VI', grid(P, T), BASE, {'P': P, 'T': T})

    # ---- Ice VII_X_French (LBF path in MATLAB via OpenGval) ----
    P = np.array([3000.0, 5000.0, 10000.0])
    T = np.array([400.0, 500.0, 700.0])
    add_case(cases, 'VII_X_French_grid', 'VII_X_French', grid(P, T), BASE,
             {'P': P, 'T': T})

    # ---- Liquid water ----
    P = np.arange(0.1, 2000.1, 500.0)
    T = np.arange(260.0, 460.0, 50.0)
    add_case(cases, 'water1_grid', 'water1', grid(P, T), BASE, {'P': P, 'T': T})
    Ps = np.array([100.0, 500.0, 1500.0])
    Ts = np.array([280.0, 320.0, 380.0])
    add_case(cases, 'water1_scatter', 'water1', scatter(Ps, Ts), BASE,
             {'P': Ps, 'T': Ts})

    P = np.array([1.0, 50.0, 200.0])
    T = np.array([280.0, 320.0, 360.0])
    add_case(cases, 'water_IAPWS95_grid', 'water_IAPWS95', grid(P, T), BASE,
             {'P': P, 'T': T})

    # water2 — Brown 2018 EOS, valid up to 100 GPa.
    P = np.array([100.0, 1000.0, 5000.0])
    T = np.array([300.0, 500.0, 800.0])
    add_case(cases, 'water2_grid', 'water2', grid(P, T), BASE, {'P': P, 'T': T})

    # ---- water3: Helmholtz fluid (vapour, liquid, supercritical) ----
    P = np.array([0.1, 10.0, 100.0, 1000.0])
    T = np.array([300.0, 400.0, 500.0, 700.0])
    add_case(cases, 'water3_grid', 'water3', grid(P, T), BASE, {'P': P, 'T': T})
    Ps = np.array([1e-3, 0.101325, 30.0, 500.0, 2000.0])        # vapour, liquid, supercritical, ...
    Ts = np.array([300.0, 298.15, 700.0, 350.0, 450.0])
    add_case(cases, 'water3_scatter', 'water3', scatter(Ps, Ts), BASE, {'P': Ps, 'T': Ts})

    # ---- (rho,T) input: Helmholtz water3 (direct) and Gibbs water1 (via rho2P) ----
    R = np.array([1e-2, 997.0, 1100.0, 1250.0])
    Tr = np.array([400.0, 298.15, 300.0, 350.0])
    add_case(cases, 'water3rhoT_scatter', 'water3', scatter(R, Tr), BASE, {'P': R, 'T': Tr}, rhoT=True)
    R = np.array([1000.0, 1050.0, 1100.0])
    Tr = np.array([280.0, 300.0, 330.0])
    add_case(cases, 'water1rhoT_grid', 'water1', grid(R, Tr), BASE, {'P': R, 'T': Tr}, rhoT=True)

    # ---- Aqueous NaCl (3D) ----
    P = np.arange(0.1, 500.1, 100.0)
    T = np.arange(273.0, 401.0, 25.0)
    m = np.array([0.1, 0.5, 1.0, 3.0])
    add_case(cases, 'NaClaq_grid', 'NaClaq', grid(P, T, m), NACL,
             {'P': P, 'T': T, 'm': m})

    # NaClaq edges: m at the cutoff (2e-4), middle, and near upper knot (~7).
    P = np.array([1.0, 200.0, 500.0])
    T = np.array([275.0, 300.0, 380.0])
    m = np.array([2e-4, 0.5, 5.0])
    add_case(cases, 'NaClaq_edges_grid', 'NaClaq', grid(P, T, m), NACL,
             {'P': P, 'T': T, 'm': m})

    # NaClaq scatter — locks in the recently fixed scatter mixing path.
    Ps = np.array([100.0, 200.0, 500.0])
    Ts = np.array([280.0, 300.0, 350.0])
    ms = np.array([0.5, 1.5, 3.0])
    add_case(cases, 'NaClaq_scatter', 'NaClaq', scatter(Ps, Ts, ms), NACL,
             {'P': Ps, 'T': Ts, 'm': ms})

    out_path = os.path.join(HERE, 'reference_SeaFreeze.mat')
    savemat(out_path, cases, do_compression=True)
    print(f"\nWrote {len(cases)} cases to {out_path}")


if __name__ == '__main__':
    main()
