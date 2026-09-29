#!/usr/bin/env python3
"""Convert SeaFreeze's MATLAB v7.3 spline files into Octave-readable MAT v7.

Ten of the spline files under ``Matlab/splines/`` were saved as MATLAB v7.3,
which is HDF5 underneath.  In that format a MATLAB cell array (``sp.knots``,
``sp.PTmc``) is stored as a dataset of HDF5 object references pointing into a
hidden ``#refs#`` group.  GNU Octave's ``load`` cannot follow those references,
so the splines come back unusable.

This script reads each v7.3 file with h5py, resolves the references itself, and
re-saves the ``sp`` struct as compressed MAT v7 (a v5-format file), which both
Octave and MATLAB read natively.  Output goes to ``Matlab/splines_octave/``,
mirroring the ``splines/`` layout; the original files are never modified.

Only the ten v7.3 files are converted.  The five NaCl_aq files are already v5
and load fine in Octave, so ``sf_load_spline`` falls back to ``splines/`` for
them rather than carrying a pointless duplicate.

Usage:
    python3 tools/convert_splines_for_octave.py            # convert
    python3 tools/convert_splines_for_octave.py --check    # verify, write nothing

Requires: h5py, scipy, numpy.
"""

import argparse
import os
import sys

import numpy as np

try:
    import h5py
    import scipy.io as sio
except ImportError as exc:  # pragma: no cover - environment problem, not logic
    sys.exit("error: this script needs h5py and scipy installed (%s)" % exc)


# (subfolder, filename) for every spline that sf_load_spline can load and that
# is stored as v7.3.  Kept in the same order as sf_load_spline's lookup table.
V73_SPLINES = [
    ("ice_Ih", "ice_Ih.mat"),
    ("ice_II", "ice_II.mat"),
    ("ice_III", "ice_III.mat"),
    ("ice_V", "ice_V.mat"),
    ("ice_VI", "ice_VI.mat"),
    ("ice_VII_X_French", "ice_VII_X_French.mat"),
    ("water_Bollengier", "water_Bollengier.mat"),
    ("water_Brown", "water_Brown.mat"),
    ("water_IAPWS95", "water_IAPWS95.mat"),
    ("NaCl_aq_Brown2024", "NaCl_aq_Brown2024.mat"),
]


def _field_names(group):
    """Decode a struct group's MATLAB_fields attribute into a list of names."""
    names = []
    for entry in group.attrs["MATLAB_fields"]:
        if isinstance(entry, np.ndarray):
            names.append(b"".join(entry).decode())
        else:
            names.append(str(entry))
    return names


def _convert(h5, obj, skipped, path):
    """Recursively turn an h5py node from a v7.3 MAT file into numpy/py data.

    Returns None for anything that cannot be represented without MATLAB's
    object subsystem (datetime and other MCOS classes); the caller drops those
    fields and records them in `skipped`.
    """
    # MCOS-backed classes (datetime, string, table, ...) live in #subsystem#
    # and carry MATLAB_object_decode.  We cannot reconstruct them here.
    if "MATLAB_object_decode" in obj.attrs:
        skipped.append("%s (%s)" % (path, obj.attrs.get("MATLAB_class", b"?").decode()))
        return None

    if isinstance(obj, h5py.Group):
        if "MATLAB_fields" not in obj.attrs:
            # A plain group with no field list: treat its members as fields.
            return {k: _convert(h5, v, skipped, "%s.%s" % (path, k))
                    for k, v in obj.items()}
        out = {}
        for name in _field_names(obj):
            value = _convert(h5, obj[name], skipped, "%s.%s" % (path, name))
            if value is not None:
                out[name] = value
        return out

    mat_class = obj.attrs.get("MATLAB_class", b"").decode()

    # Empty arrays / structs are flagged rather than stored with a real shape.
    if obj.attrs.get("MATLAB_empty", 0):
        return "" if mat_class == "char" else np.zeros((0, 0))

    if mat_class == "cell":
        raw = obj[()]
        # HDF5 stores MATLAB arrays with reversed axes, cells included.
        out = np.empty(raw.shape[::-1], dtype=object)
        for idx in np.ndindex(raw.shape):
            out[idx[::-1]] = _convert(h5, h5[raw[idx]], skipped,
                                      "%s{%s}" % (path, ",".join(map(str, idx))))
        return out

    if mat_class == "char":
        return "".join(chr(c) for c in np.asarray(obj[()]).ravel())

    arr = np.asarray(obj[()])
    # Undo the HDF5 axis reversal so coefs comes back as n1 x n2 x ... and not
    # its transpose.  1-D data needs no reordering.
    return arr.T if arr.ndim > 1 else arr


def _validate(sp, label):
    """Check the B-form invariants that sp_val relies on.  Returns a message list."""
    problems = []
    for key in ("knots", "coefs", "number", "order"):
        if key not in sp:
            problems.append("missing field '%s'" % key)
    if problems:
        return problems

    number = np.asarray(sp["number"]).ravel().astype(int)
    order = np.asarray(sp["order"]).ravel().astype(int)
    knots = np.asarray(sp["knots"]).ravel()
    coefs = np.asarray(sp["coefs"])

    if len(knots) != len(number):
        problems.append("%d knot vectors but number has %d entries"
                        % (len(knots), len(number)))
    else:
        for i, kv in enumerate(knots):
            expected = number[i] + order[i]
            if np.size(kv) != expected:
                problems.append("dim %d: len(knots)=%d, expected number+order=%d"
                                % (i + 1, np.size(kv), expected))

    if coefs.shape != tuple(number):
        problems.append("coefs shape %s does not match number %s"
                        % (coefs.shape, tuple(number)))

    if not problems:
        print("    ok: order=%s number=%s coefs=%s knots=%s"
              % ([int(v) for v in order], [int(v) for v in number],
                 tuple(int(v) for v in coefs.shape),
                 [int(np.size(k)) for k in knots]))
    return problems


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--check", action="store_true",
                        help="validate only; do not write any files")
    args = parser.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    matlab_dir = os.path.dirname(here)
    src_root = os.path.join(matlab_dir, "splines")
    dst_root = os.path.join(matlab_dir, "splines_octave")

    failures = 0
    for subfolder, filename in V73_SPLINES:
        src = os.path.join(src_root, subfolder, filename)
        dst = os.path.join(dst_root, subfolder, filename)
        print("%s/%s" % (subfolder, filename))

        if not os.path.isfile(src):
            print("    ERROR: source not found: %s" % src)
            failures += 1
            continue

        skipped = []
        with h5py.File(src, "r") as h5:
            if "sp" not in h5:
                print("    ERROR: no 'sp' variable (not a SeaFreeze spline?)")
                failures += 1
                continue
            sp = _convert(h5, h5["sp"], skipped, "sp")

        for problem in _validate(sp, filename):
            print("    ERROR: %s" % problem)
            failures += 1
        if skipped:
            # datetime provenance fields and the like: metadata fnGval never reads.
            print("    dropped (MATLAB object, unreadable outside MATLAB): %s"
                  % ", ".join(skipped))

        if args.check:
            continue

        os.makedirs(os.path.dirname(dst), exist_ok=True)
        sio.savemat(dst, {"sp": sp}, format="5", do_compression=True)
        print("    wrote %s (%.1f kB)" % (os.path.relpath(dst, matlab_dir),
                                          os.path.getsize(dst) / 1024.0))

    print()
    if failures:
        print("%d problem(s) found" % failures)
        return 1
    print("%d spline(s) %s" % (len(V73_SPLINES),
                               "checked" if args.check else "converted"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
