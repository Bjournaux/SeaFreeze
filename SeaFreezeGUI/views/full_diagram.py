"""Full H2O phase diagram (water3 + ices Ih-VI) in P-T or rho-T, coloured by
stable phase or by any property of the stable phase (2-D map, map +
isocontours, or 3-D surface), with the boundaries, saturation curve,
critical point and triple points overlaid."""

import numpy as np
import pandas as pd
import streamlit as st

from core import diagram_plot as P
from core import diagrams as D
from core.constants import PROP_PRESSURE, PROPS_ALL, PROPS_SOLID
from core.ui import make_csv

_COORDS = {"P–T": "PT", "ρ–T": "rhoT"}
_VIEWS = ["2-D map", "2-D map + isocontours", "3-D surface"]
_PHASE = "Stable phase"
_PROP_META = {**PROP_PRESSURE, **PROPS_ALL, **PROPS_SOLID}
_COLORSCALES = ["Viridis", "Cividis", "Plasma", "Inferno", "Magma", "Turbo", "YlGnBu", "RdBu",
                "Spectral"]
_OVERLAY_LABELS = {"boundaries": "Phase boundaries", "saturation": "Saturation curve / dome",
                   "critical": "Critical point", "triple": "Triple points / tie lines",
                   "labels": "Field labels"}


def _prop_label(p):
    name, unit = _PROP_META[p]
    return f"{p} — {name}" + (f" ({unit})" if unit != "—" else "")


def _applied(coords):
    """The (window, resolution) currently applied for these coordinates."""
    store = st.session_state.setdefault("pdf_applied", {})
    if coords not in store:
        store[coords] = (dict(D.DEFAULT_WINDOW[coords]), D.DEFAULT_RESOLUTION)
    return store[coords]


def _range_form(coords):
    """Sidebar form for the window and resolution; applies on submit only."""
    window, res = _applied(coords)
    default = D.is_default(coords, window, res)
    with st.sidebar.expander("Range & resolution", expanded=not default):
        with st.form(f"pdf_form_{coords}"):
            if coords == "PT":
                c1, c2 = st.columns(2)
                p0 = c1.number_input("P min (MPa)", value=float(window["P"][0]), min_value=1e-12,
                                     max_value=1e7, format="%.3g")
                p1 = c2.number_input("P max (MPa)", value=float(window["P"][1]), min_value=1e-12,
                                     max_value=1e7, format="%.3g")
            else:
                c1, c2 = st.columns(2)
                r0 = c1.number_input("ρ min (kg/m³)", value=float(window["rho"][0]), min_value=1e-12,
                                     max_value=16000.0, format="%.3g")
                r1 = c2.number_input("ρ max (kg/m³)", value=float(window["rho"][1]), min_value=1e-12,
                                     max_value=16000.0, format="%.3g")
                xs = st.radio("Density grid spacing", ["log", "linear"], horizontal=True,
                              index=0 if window["xscale"] == "log" else 1)
            c1, c2 = st.columns(2)
            t0 = c1.number_input("T min (K)", value=float(window["T"][0]), min_value=1.0,
                                 max_value=1e5, step=10.0)
            t1 = c2.number_input("T max (K)", value=float(window["T"][1]), min_value=1.0,
                                 max_value=1e5, step=10.0)
            new_res = st.select_slider("Resolution", list(D.RESOLUTIONS), value=res)
            st.caption("⚠ Any range or resolution other than the default is computed live: "
                       + ", ".join(f"{k} ≈ {v[coords]} s" for k, v in D.RUNTIME_S.items())
                       + " on a recent laptop — longer on slower machines. "
                         "Results are cached for the session.")
            go_btn = st.form_submit_button("Recompute", type="primary", use_container_width=True)
        if st.button("Reset to default (instant)", use_container_width=True,
                     key=f"pdf_reset_{coords}", disabled=default):
            st.session_state["pdf_applied"][coords] = (dict(D.DEFAULT_WINDOW[coords]),
                                                       D.DEFAULT_RESOLUTION)
            st.rerun()
    if go_btn:
        if coords == "PT":
            win = dict(P=(min(p0, p1), max(p0, p1)), T=(min(t0, t1), max(t0, t1)))
        else:
            # the underlying P-T map spans the default pressures
            win = dict(rho=(min(r0, r1), max(r0, r1)), T=(min(t0, t1), max(t0, t1)),
                       P=D.DEFAULT_WINDOW["rhoT"]["P"], xscale=xs)
        bad = (win["T"][0] == win["T"][1]
               or (coords == "PT" and win["P"][0] == win["P"][1])
               or (coords == "rhoT" and win["rho"][0] == win["rho"][1]))
        if bad:
            st.sidebar.error("Each range needs two different end values.")
        else:
            st.session_state["pdf_applied"][coords] = (win, new_res)
            st.rerun()


def render():
    with st.sidebar:
        st.header("Full diagram")
        coords = _COORDS[st.radio("Coordinates", list(_COORDS), horizontal=True, key="pdf_coords")]
        props = [_PHASE] + list(D.MAP_PROPS) + list(D.ICE_ONLY_PROPS)
        colour = st.selectbox("Colour by", props, key="pdf_colour",
                              format_func=lambda p: p if p == _PHASE else _prop_label(p))
        is_prop = colour != _PHASE
        view = st.radio("View", _VIEWS if is_prop else _VIEWS[:1], key="pdf_view",
                        help=None if is_prop else "Isocontours and the 3-D surface need a property.")
        colorscale, logz, ncont, brk, crange = "Viridis", False, 15, True, "1–99 %"
        if is_prop:
            c1, c2 = st.columns(2)
            colorscale = c1.selectbox("Colour scale", _COLORSCALES, key="pdf_cscale")
            logz = c2.checkbox("Log colour scale", value=colour in P.LOG_DEFAULT,
                               key=f"pdf_log_{colour}",
                               help="log10 of the property (positive values only).")
            crange = st.radio("Colour range", ["1–99 %", "Full", "Manual"], horizontal=True,
                              key="pdf_crange",
                              help="1–99 %: the extreme 1 % at each end (e.g. Cp at the critical "
                                   "point) is shown in the end colours, so the rest of the map "
                                   "keeps its contrast.")
            if crange == "Manual":
                c1, c2 = st.columns(2)
                u = _PROP_META[colour][1]
                cmin = c1.number_input(f"min ({u})", value=0.0, format="%.4g", key=f"pdf_cmin_{colour}")
                cmax = c2.number_input(f"max ({u})", value=1.0, format="%.4g", key=f"pdf_cmax_{colour}")
            if view == _VIEWS[1]:
                ncont = st.slider("Isocontours", 5, 30, 15, key="pdf_ncont")
            if view == _VIEWS[2]:
                brk = st.checkbox("Open the surface at phase boundaries", value=True, key="pdf_break",
                                  help="Leaves a gap along each boundary instead of a vertical wall.")
        st.subheader("Axes")
        c1, c2 = st.columns(2)
        xlog = c1.checkbox("log P" if coords == "PT" else "log ρ", key=f"pdf_xlog_{coords}",
                           value=coords == "PT" or _applied(coords)[0]["xscale"] == "log")
        ylog = c2.checkbox("log T", value=False, key="pdf_ylog")
        st.subheader("Overlays")
        overlays = [k for k, lab in _OVERLAY_LABELS.items()
                    if st.checkbox(lab, value=True, key=f"pdf_ov_{k}")]
    _range_form(coords)

    window, res = _applied(coords)
    key = D.window_items(window)
    default = D.is_default(coords, window, res)
    if default:
        diag = D.get_diagram(coords, key, res)
    else:
        with st.spinner(f"Computing the {'P–T' if coords == 'PT' else 'ρ–T'} diagram at {res} "
                        f"resolution (≈ {D.RUNTIME_S[res][coords]} s)…"):
            diag = D.get_diagram(coords, key, res)
    pmap = None
    if is_prop:
        with st.spinner("Evaluating the property maps (a few seconds, once per diagram)…"):
            pmap = D.get_property_maps(coords, key, res)

    # ── caption ───────────────────────────────────────────────────────────
    xs = ("P " + " – ".join(f"{v:.3g}" for v in window["P"]) + " MPa") if coords == "PT" else \
         ("ρ " + " – ".join(f"{v:.3g}" for v in window["rho"]) + " kg/m³")
    src = "precomputed default" if default else "computed live"
    st.caption(f"Stable phase by Gibbs-energy minimisation — fluid **water3** (Helmholtz, beta), "
               f"ices **Ih, II, III, V, VI**. {xs}, T {window['T'][0]:.0f} – {window['T'][1]:.0f} K, "
               f"{res} resolution ({src})."
               + (" Grey: two-phase regions." if coords == "rhoT" else ""))

    unit = _PROP_META[colour][1] if is_prop else ""
    cr = {"1–99 %": "robust", "Full": "full"}.get(crange)
    if cr is None:                                       # manual, in data units
        lo, hi = sorted((cmin, cmax))
        if logz and lo <= 0:
            st.warning("Log colour scale: the manual minimum must be positive; using 1–99 %.")
            cr = "robust"
        else:
            cr = (np.log10(lo), np.log10(hi)) if logz else (lo, hi)
    if is_prop and not np.isfinite(pmap.values[colour]).any():
        st.warning(f"{colour} is not defined anywhere in this window.")
        return
    if view == _VIEWS[2]:
        fig = P.figure_3d(diag, pmap, colour, unit, overlays, colorscale, logz, brk,
                          xlog=xlog, ylog=ylog, crange=cr)
    else:
        fig = P.figure_2d(diag, pmap, colour if is_prop else None, unit, overlays, colorscale,
                          logz, contours=view == _VIEWS[1], ncontours=ncont, xlog=xlog, ylog=ylog,
                          crange=cr)
    st.plotly_chart(fig, use_container_width=True, key=f"pdf_chart_{coords}")

    notes = ["Ice VII/X is not included yet (a new representation is in preparation): the "
             "high-pressure region above ice VI is left blank."
             if (diag.pm.stable < 0).any() else None,
             "The fluid is not used more than 40 K below the melting curve of the stable solid "
             "(the surface's validity mask); below 230 K the vapour is the surface's ideal-gas "
             "part (dilute-vapour extension)."]
    if is_prop and colour in D.ICE_ONLY_PROPS:
        notes.append(f"{colour} is defined for the ices only (blank in the fluid).")
    for n in filter(None, notes):
        st.caption("• " + n)

    _details(diag, pmap, colour if is_prop else None, coords)


def _details(diag, pmap, prop, coords):
    """Triple points / critical point table and CSV export of the displayed grid."""
    tps = diag.triple_points if coords == "PT" else diag.tie_lines
    with st.expander("Triple points and critical point"):
        rows = [dict(Point=" – ".join(t["labels"]), **{"P (MPa)": t["P"], "T (K)": t["T"]},
                     **{f"ρ {l} (kg/m³)": r for l, r in zip(("1", "2", "3"), t["rho"])})
                for t in tps]
        c = diag.critical
        rows.append({"Point": "critical point", "P (MPa)": c["P"], "T (K)": c["T"],
                     "ρ 1 (kg/m³)": c["rho"]})
        st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
        st.caption("Triple points refined by Newton on G_a = G_b = G_c; ρ 1–3 are the densities "
                   "of the three coexisting phases, in the order of the point's name.")

    x = diag.pm.P if coords == "PT" else diag.rho
    field = diag.pm.stable if coords == "PT" else diag.field
    X, T = np.meshgrid(x, diag.pm.T if coords == "PT" else diag.T, indexing="ij")
    names = np.array(list(diag.pm.names) + ["two-phase"], dtype=object)
    df = pd.DataFrame({"P (MPa)" if coords == "PT" else "rho (kg/m3)": X.ravel(), "T (K)": T.ravel(),
                       "phase": np.where(field.ravel() < 0, "", names[np.clip(field.ravel(), 0, None)])})
    if prop is not None:
        df[f"{prop} ({_PROP_META[prop][1]})"] = pmap.values[prop].ravel()
    st.download_button("Download grid (CSV)", make_csv(df),
                       file_name=f"seafreeze_phase_diagram_{coords}{'_' + prop if prop else ''}.csv",
                       mime="text/csv", key=f"pdf_csv_{coords}")
