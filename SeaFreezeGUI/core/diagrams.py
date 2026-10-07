"""Full H2O phase diagram data for the GUI: precomputed defaults + cached live runs.

The diagrams come from seafreeze.phasediagram (phase_diagram_PT /
phase_diagram_rhoT, fluid water_Brown2026 + ices Ih-VI) and the property maps from
seafreeze.phasediagram.property_map.  The default window at the default
resolution is shipped precomputed in assets/diagrams/ (regenerate with
tools/precompute_diagrams.py) so the page opens instantly; any other window
or resolution is computed live and cached for the session.
"""

import json
import os
import warnings

import numpy as np
import streamlit as st

from seafreeze import __version__ as _sf_version
from seafreeze import phasediagram as _pd
from seafreeze.coexistence import Coexistence

from core.ui import asset

FLUID = "water_Brown2026"
ICES = tuple(_pd.ICES)

DEFAULT_WINDOW = {
    "PT":   dict(P=(1e-8, 1e5), T=(150.0, 1800.0)),
    "rhoT": dict(rho=(1e-7, 4e3), T=(150.0, 1800.0), P=(1e-10, 1e5), xscale="log"),
}
RESOLUTIONS = {
    "Draft":    dict(PT=dict(nP=200, nT=160), rhoT=dict(nP=500, nT=160, nrho=320)),
    "Standard": dict(PT=dict(nP=400, nT=320), rhoT=dict(nP=1000, nT=300, nrho=600)),
    "Fine":     dict(PT=dict(nP=700, nT=560), rhoT=dict(nP=1600, nT=480, nrho=960)),
}
DEFAULT_RESOLUTION = "Standard"
# Typical live runtime (s) on an Apple M-series laptop, diagram + all property maps
RUNTIME_S = {
    "Draft":    dict(PT=8, rhoT=9),
    "Standard": dict(PT=12, rhoT=15),
    "Fine":     dict(PT=23, rhoT=26),
}

# bump when the diagram data change without a SeaFreeze version change, so
# stale precomputed files are recomputed rather than loaded
DATA_VERSION = 5

MAP_PROPS = tuple(_pd.MAP_PROPS)
ICE_ONLY_PROPS = tuple(_pd.ICE_ONLY_PROPS)


def window_items(window):
    """Hashable, canonical form of a window dict (cache key)."""
    return tuple(sorted((k, tuple(float(x) for x in v) if isinstance(v, (tuple, list)) else v)
                        for k, v in window.items()))


def is_default(coords, window, resolution):
    return (resolution == DEFAULT_RESOLUTION
            and window_items(window) == window_items(DEFAULT_WINDOW[coords]))


def precomputed_path(coords):
    return asset(os.path.join("diagrams", f"phase_diagram_{coords}.npz"))


# ─────────────────────────────────────────────────────────────────────────────
# Build / load
# ─────────────────────────────────────────────────────────────────────────────
def build_diagram(coords, window, resolution):
    """Compute a diagram (no caching)."""
    grid = RESOLUTIONS[resolution][coords]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if coords == "PT":
            d = _pd.phase_diagram_PT(P=window["P"], T=window["T"], fluid=FLUID, ices=ICES, **grid)
        else:
            d = _pd.phase_diagram_rhoT(rho=window["rho"], T=window["T"], P=window["P"],
                                       xscale=window["xscale"], fluid=FLUID, ices=ICES, **grid)
    d.pm.G = None       # every phase's G on the grid: not used by the maps; frees memory
    return d


# The caches are shared by every session of a deployed app (~1 GB on Streamlit
# Community Cloud): keep a few diagrams only.
@st.cache_resource(show_spinner=False, max_entries=4)
def get_diagram(coords, window_key, resolution):
    """The diagram for a window (window_items form) and resolution.

    The shipped default loads from disk; anything else is computed.  Cached
    as a shared resource (these are large; property_map does not mutate them).
    """
    window = dict(window_key)
    if is_default(coords, window, resolution):
        d = load_diagram(precomputed_path(coords), coords, window, resolution)
        if d is not None:
            return d
    return build_diagram(coords, window, resolution)


@st.cache_resource(show_spinner=False, max_entries=4)
def get_property_maps(coords, window_key, resolution):
    """PropertyMap with every MAP_PROPS + ICE_ONLY_PROPS map of the diagram."""
    d = get_diagram(coords, window_key, resolution)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return _pd.property_map(d, *(MAP_PROPS + ICE_ONLY_PROPS))


# ─────────────────────────────────────────────────────────────────────────────
# (De)serialisation of the precomputed defaults — plain arrays + JSON, no pickle
# ─────────────────────────────────────────────────────────────────────────────
def _meta(coords, window, resolution):
    return dict(coords=coords, window=window_items(window), resolution=resolution,
                grid=RESOLUTIONS[resolution][coords], fluid=FLUID, ices=list(ICES),
                seafreeze=_sf_version, data_version=DATA_VERSION)


def _jsonable(x):
    if isinstance(x, dict):
        return {k: _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (np.floating, float)):
        return float(x)
    if isinstance(x, (np.integer, int)):
        return int(x)
    return x


def _pack_pm(pm, out):
    out["pm_P"], out["pm_T"] = pm.P, pm.T
    out["pm_stable"] = pm.stable.astype(np.int16)
    # each phase's density only where it is stable (all that the maps use)
    rho = np.where(pm.stable[None] == np.arange(len(pm.names))[:, None, None], pm.rho, np.nan)
    out["pm_rho"] = rho.astype(np.float32)
    out["pm_rho_stable"] = pm.rho_stable.astype(np.float32)


def _unpack_pm(z, names):
    return _pd.PhaseMap(P=z["pm_P"], T=z["pm_T"], names=names, G=None,
                        rho=z["pm_rho"].astype(float), stable=z["pm_stable"].astype(int),
                        rho_stable=z["pm_rho_stable"].astype(float))


def _pack_lines(prefix, lines, out):
    """Ragged list of equal-length array tuples -> index + concatenated arrays."""
    idx = [(int(l[0]), int(l[1]), len(l[2])) for l in lines]
    out[prefix + "_idx"] = np.array(idx, dtype=np.int64).reshape(-1, 3)
    for k in range(2, len(lines[0]) if lines else 2):
        out[f"{prefix}_{k}"] = np.concatenate([np.asarray(l[k], float) for l in lines])


def _unpack_lines(prefix, z, width):
    idx = z[prefix + "_idx"]
    cols = [z[f"{prefix}_{k}"] for k in range(2, width)] if len(idx) else []
    out, a = [], 0
    for i, j, n in idx:
        out.append((int(i), int(j), *[c[a:a + n] for c in cols]))
        a += n
    return out


def _pack_sat(s, out):
    if s is not None:
        for f in Coexistence._fields:
            out["sat_" + f] = np.asarray(getattr(s, f), float)


def _unpack_sat(z):
    if "sat_P" not in z:
        return None
    return Coexistence(*[z["sat_" + f] for f in Coexistence._fields])


def save_diagram(d, path, coords, window, resolution):
    out = {}
    _pack_pm(d.pm, out)
    _pack_sat(d.saturation, out)
    info = dict(meta=_meta(coords, window, resolution), names=list(d.pm.names),
                critical=_jsonable(d.critical))
    if coords == "PT":
        _pack_lines("bnd", d.boundaries, out)
        info.update(triple_points=_jsonable(d.triple_points), labels=_jsonable(d.labels))
    else:
        out["rho"], out["T"] = d.rho, d.T
        out["field"] = d.field.astype(np.int16)
        out["gap"] = d.gap.astype(np.int8)
        _pack_lines("coex", d.coexistence, out)
        info.update(two_phase=int(d.two_phase), tie_lines=_jsonable(d.tie_lines),
                    labels=_jsonable(d.labels), xscale=d.xscale)
    out["info"] = np.array(json.dumps(info))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez_compressed(path, **out)


def load_diagram(path, coords, window, resolution):
    """The saved diagram if it matches (coords, window, resolution, library), else None."""
    if not os.path.isfile(path):
        return None
    try:
        with np.load(path, allow_pickle=False) as f:
            z = {k: f[k] for k in f.files}
        info = json.loads(str(z["info"]))
        want = json.loads(json.dumps(_meta(coords, window, resolution)))
        if info["meta"] != want:
            return None
        names = info["names"]
        pm = _unpack_pm(z, names)
        sat = _unpack_sat(z)
        tp = [dict(t, rho=tuple(t["rho"])) for t in info.get("triple_points", info.get("tie_lines", []))]
        labels = [tuple(l) for l in info["labels"]]
        if coords == "PT":
            return _pd.DiagramPT(pm=pm, boundaries=_unpack_lines("bnd", z, 4), saturation=sat,
                                 critical=info["critical"], triple_points=tp, labels=labels,
                                 missing=pm.stable < 0)
        return _pd.DiagramRhoT(pm=pm, rho=z["rho"], T=z["T"], field=z["field"].astype(int),
                               two_phase=info["two_phase"], coexistence=_unpack_lines("coex", z, 5),
                               saturation=sat, critical=info["critical"], tie_lines=tp,
                               labels=labels, xscale=info["xscale"],
                               gap=z["gap"].astype(np.int8) if "gap" in z else None)
    except Exception:
        return None                                  # stale or unreadable: recompute
