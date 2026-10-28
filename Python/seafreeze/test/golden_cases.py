"""Shared case definitions for the golden-value regression suite.

Each material gets three physically meaningful (P, T[, m]) axes.  The
scatter case zips the axes point-wise (3 points); the grid case takes their
outer product (3x3 or 3x3x3).  Both are evaluated in all-properties mode.

Edit this file to add materials or points, then regenerate the reference with

    python -m seafreeze.test.gen_golden_reference

and commit the updated ``golden_reference.json`` together with the change
that motivated it.
"""

# material -> dict(P=[...], T=[...], m=[...] or None)
GOLDEN_AXES = {
    'Ih':             dict(P=[0.1,   100.0,   200.0],   T=[200.0, 240.0,  265.0],  m=None),
    'II':             dict(P=[200.0, 400.0,   600.0],   T=[180.0, 220.0,  250.0],  m=None),
    'III':            dict(P=[220.0, 300.0,   400.0],   T=[230.0, 245.0,  255.0],  m=None),
    'V':              dict(P=[400.0, 600.0,   800.0],   T=[230.0, 250.0,  270.0],  m=None),
    'VI':             dict(P=[800.0, 1500.0,  2500.0],  T=[230.0, 280.0,  350.0],  m=None),
    'VII_X_French':   dict(P=[2500.0, 10000.0, 50000.0], T=[300.0, 600.0, 1200.0], m=None),
    'water1':         dict(P=[0.1,   500.0,   1500.0],  T=[260.0, 300.0,  400.0],  m=None),
    # Dense supercritical fluid only: Brown 2018 is unphysical in the low-T
    # corners where ices are stable (negative Cp, non-monotonic rho).
    'water2':         dict(P=[100.0, 2000.0,  20000.0], T=[600.0, 1500.0, 4000.0], m=None),
    'water_IAPWS95':  dict(P=[0.1,   100.0,   1000.0],  T=[280.0, 400.0,  800.0],  m=None),
    # Stitched: one point LP-only, one in the taper, one HP-only.
    'NaClaq':         dict(P=[100.0, 750.0,   3000.0],  T=[250.0, 300.0,  450.0],  m=[0.5, 2.0, 5.0]),
    'NaClaq_LP':      dict(P=[10.0,  300.0,   900.0],   T=[250.0, 300.0,  450.0],  m=[0.5, 2.0, 5.0]),
    'NaClaq_HP':      dict(P=[600.0, 2000.0,  8000.0],  T=[300.0, 600.0,  1500.0], m=[0.5, 2.0, 5.0]),
    'NaClaq_5GPa_2024': dict(P=[100.0, 1000.0, 4000.0], T=[250.0, 300.0,  450.0],  m=[0.5, 2.0, 5.0]),
}

# Relative tolerance for the comparison.  Spline evaluation is deterministic;
# differences across NumPy/SciPy versions are at the 1e-13 level.  1e-8 leaves
# room for BLAS reordering while still catching any real change in a formula,
# a constant, or a spline coefficient.
GOLDEN_RTOL = 1e-8
