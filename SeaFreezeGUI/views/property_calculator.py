"""Property Calculator page — single points, 1-D sweeps and 2-D grids.

Input in (P, T[, m]) for every material, or (rho, T) for the pure phases
(native for the Helmholtz fluid water_Brown2026; by pressure inversion for the Gibbs
phases).  Every P, rho and T range can be log-spaced and drawn on a log axis.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from core.constants import (
    ALL_MATERIALS, INPUT_LIMITS, available_properties, categorized_properties, is_helmholtz,
    is_nacl, supports_rhoT,
)
from core.compute import (
    compute_properties, get_phase_range, get_rho_range, get_stability_boundaries,
)
from core.ui import (
    BOUNDARY_COLORS, add_phase_boundaries, fmt_range, label, make_csv, title_with_fixed,
)

# display name of each coordinate (plot axes, titles, CSV headers)
_VAR_LABEL = {"P": "P", "rho": "ρ", "T": "T", "m": "m"}
_BRANCHES = {"stable": "stable (lower G)", "liquid": "liquid", "vapor": "vapour"}
# (rho, T) input for a Gibbs phase inverts P at each point
_GIBBS_RHOT_MS = 8                  # ~ms per point (Apple M-series)
_GIBBS_RHOT_MAX = 2500              # points


def _axis_title(name, unit):
    return f"{_VAR_LABEL.get(name, name)} ({unit})"


def _fmt(lo, hi, unit):
    """Range text that stays readable from 1e-10 to 1e7."""
    if lo > 0 and (hi / lo > 1e4 or hi >= 1e5 or lo < 0.1):
        return f"{lo:.3g} – {hi:.3g} {unit}"
    return fmt_range(lo, hi, unit)


def _coord_input(var, title, unit, bounds, rng, single, step, log_default, ctx):
    """Single value or range for one coordinate.

    Returns (values, is_range, log).  Value widgets are keyed per material and
    input mode (ctx) so each material starts from its own defaults.
    """
    lo, hi = (float(b) for b in bounds)
    st.subheader(f"{title} ({unit})")
    mode = st.radio(f"{var} input", ["Single value", "Range"], key=f"{var}_mode",
                    horizontal=True, label_visibility="collapsed")
    fmt = "%.6g"
    if mode == "Single value":
        v = st.number_input(f"{_VAR_LABEL[var]} ({unit})", value=float(np.clip(single, lo, hi)),
                            min_value=lo, max_value=hi, step=step, format=fmt,
                            key=f"{var}_single_{ctx}")
        return [v], False, False
    c1, c2 = st.columns(2)
    vmin = c1.number_input(f"{_VAR_LABEL[var]} min", value=float(np.clip(rng[0], lo, hi)),
                           min_value=lo, max_value=hi, step=step, format=fmt,
                           key=f"{var}_min_{ctx}")
    vmax = c2.number_input(f"{_VAR_LABEL[var]} max", value=float(np.clip(rng[1], lo, hi)),
                           min_value=lo, max_value=hi, step=step, format=fmt,
                           key=f"{var}_max_{ctx}")
    c1, c2 = st.columns([3, 2])
    n = c1.number_input(f"{_VAR_LABEL[var]} points", value=50, min_value=2, max_value=500,
                        key=f"{var}_npts")
    can_log = vmin > 0 and vmax > 0
    log = c2.checkbox("log", value=log_default and can_log, key=f"{var}_log_{ctx}",
                      disabled=not can_log,
                      help="Log-spaced points and a log axis (needs positive bounds).")
    log = bool(log and can_log)
    vals = (np.geomspace if log else np.linspace)(vmin, vmax, int(n)).tolist()
    return vals, True, log


def render():
    st.title("SeaFreeze — Property Calculator")

    with st.sidebar:
        st.header("Material")
        material = st.selectbox("Phase / Material", ALL_MATERIALS, format_func=label,
                                key="pc_material")
        nacl = is_nacl(material)
        helm = is_helmholtz(material)
        lim = INPUT_LIMITS.get(material, {})

        P_range, T_range, m_range = get_phase_range(material)
        if lim:
            P_range, T_range = lim["P"], lim["T"]

        st.divider()
        st.header("Coordinates")
        rhoT = False
        if supports_rhoT(material):
            rhoT = st.radio("Input", ["P, T", "ρ, T"], horizontal=True, key="pc_input",
                            help="(ρ, T): density and temperature in, pressure out. Native for "
                                 "water_Brown2026; for the Gibbs phases P is solved at each point "
                                 "(slower).") == "ρ, T"
        if rhoT:
            rho_range = lim.get("rho") or get_rho_range(material)
            st.caption(f"**ρ**: {_fmt(*rho_range, 'kg/m³')}  \n**T**: {_fmt(*T_range, 'K')}")
            if not helm:
                st.caption(f"⚠ (ρ, T) input for a Gibbs phase inverts P at every point "
                           f"(≈ {_GIBBS_RHOT_MS} ms per point); grids are limited to "
                           f"{_GIBBS_RHOT_MAX} points.")
        else:
            st.caption(f"**P**: {_fmt(*P_range, 'MPa')}  \n**T**: {_fmt(*T_range, 'K')}"
                       + (f"  \n**m**: {fmt_range(*m_range, 'mol/kg')}" if nacl else ""))
        ctx = f"{material}_{'rhoT' if rhoT else 'PT'}"

        # ── first coordinate: pressure or density ─────────────────────────
        if rhoT:
            rr = lim.get("rho_default", rho_range)
            X_values, X_is_range, X_log = _coord_input(
                "rho", "Density", "kg/m³", rho_range, rr,
                lim.get("rho_single", 0.5 * (rho_range[0] + rho_range[1])), 1.0,
                log_default=helm, ctx=ctx)
            X_name, X_unit = "rho", "kg/m³"
        else:
            X_values, X_is_range, X_log = _coord_input(
                "P", "Pressure", "MPa", P_range, lim.get("P_default", P_range),
                lim.get("P_single", round((P_range[0] + P_range[1]) / 2, 1)), 10.0,
                log_default=helm, ctx=ctx)
            X_name, X_unit = "P", "MPa"

        # ── temperature ───────────────────────────────────────────────────
        T_values, T_is_range, T_log = _coord_input(
            "T", "Temperature", "K", T_range, lim.get("T_default", T_range),
            lim.get("T_single", round((T_range[0] + T_range[1]) / 2, 1)), 5.0,
            log_default=False, ctx=ctx)

        # ── molality (NaClaq only) ────────────────────────────────────────
        m_values = None
        m_is_range = False
        if nacl:
            st.subheader("Molality (mol/kg)")
            m_mode = st.radio("m input", ["Single value", "Range"], key="m_mode",
                              horizontal=True, label_visibility="collapsed")
            if m_mode == "Single value":
                m_values = [st.number_input("m (mol/kg)", value=1.0,
                                            min_value=float(m_range[0]),
                                            max_value=float(m_range[1]),
                                            step=0.5, key="m_single")]
            else:
                mc1, mc2 = st.columns(2)
                m_min = mc1.number_input("m min", value=float(m_range[0]),
                                         min_value=float(m_range[0]),
                                         max_value=float(m_range[1]), step=0.5, key="m_min")
                m_max = mc2.number_input("m max", value=float(m_range[1]),
                                         min_value=float(m_range[0]),
                                         max_value=float(m_range[1]), step=0.5, key="m_max")
                m_npts = st.number_input("m points", value=20, min_value=2, max_value=200,
                                         key="m_npts")
                m_values = np.linspace(m_min, m_max, int(m_npts)).tolist()
                m_is_range = True

        # ── fluid branch (Helmholtz at P, T) ──────────────────────────────
        branch = "stable"
        if helm and not rhoT:
            branch = st.radio("Fluid branch", list(_BRANCHES), format_func=_BRANCHES.get,
                              horizontal=True, key="pc_branch",
                              help="Where liquid and vapour both exist at (P, T): the stable "
                                   "one (lower Gibbs energy) or a chosen, possibly metastable, "
                                   "branch. NaN where the chosen branch does not exist.")

        # ── property selection ────────────────────────────────────────────
        st.divider()
        st.header("Properties")
        avail = available_properties(material, rhoT)
        prop_labels = {k: f"{k}  ({v[0]}, {v[1]})" for k, v in avail.items()}
        select_all = st.checkbox("Compute all properties", value=True, key="select_all")
        if select_all:
            selected_props = list(avail.keys())
        else:
            selected_props = st.multiselect(
                "Select properties", options=list(avail.keys()), default=["rho", "Cp", "G"],
                format_func=lambda k: prop_labels[k], key="prop_select")

        st.divider()
        compute_btn = st.button("Compute", type="primary", use_container_width=True)

    # ── Determine dimensionality ──────────────────────────────────────────
    range_vars, const_vars, logs = [], [], {}
    for name, vals, is_rng, log, unit in ((X_name, X_values, X_is_range, X_log, X_unit),
                                          ("T", T_values, T_is_range, T_log, "K")):
        if is_rng:
            range_vars.append((name, vals, unit))
            logs[name] = log
        else:
            const_vars.append((_VAR_LABEL[name], vals[0], unit))
    if nacl:
        if m_is_range:
            range_vars.append(("m", m_values, "mol/kg"))
            logs["m"] = False
        else:
            const_vars.append(("m", m_values[0], "mol/kg"))
    ndim = len(range_vars)
    npoints = int(np.prod([len(v[1]) for v in range_vars])) if range_vars else 1

    mode_labels = {0: "Single point", 1: "1-D sweep (line plots)",
                   2: "2-D grid (heatmap / 3-D surface)"}
    if ndim <= 2:
        mode_str = f"**Mode**: {mode_labels[ndim]}"
        if ndim > 0:
            mode_str += (" — varying **"
                         + " & ".join(_VAR_LABEL[v[0]] + (" (log)" if logs[v[0]] else "")
                                      for v in range_vars) + "**")
        if const_vars:
            mode_str += " | " + ", ".join(f"{v[0]} = {v[1]:.6g} {v[2]}" for v in const_vars)
        if helm and not rhoT and branch != "stable":
            mode_str += f" | branch: {_BRANCHES[branch]}"
        st.info(mode_str)
    else:
        st.warning("3 independent ranges selected. Only 2 can vary at once for a "
                   "surface plot. Set one variable to a single value.")
    slow_gibbs = rhoT and not helm
    if slow_gibbs and npoints > 1:
        msg = (f"(ρ, T) input for {label(material)}: {npoints} points ≈ "
               f"{max(1, round(npoints * _GIBBS_RHOT_MS / 1000))} s (P is solved at each point).")
        (st.error if npoints > _GIBBS_RHOT_MAX else st.caption)(
            msg + (f" Reduce the grid to ≤ {_GIBBS_RHOT_MAX} points." if npoints > _GIBBS_RHOT_MAX
                   else ""))

    # ── Compute & store ───────────────────────────────────────────────────
    too_big = slow_gibbs and npoints > _GIBBS_RHOT_MAX
    if compute_btn and selected_props and ndim <= 2 and not too_big:
        try:
            mode = "scatter" if ndim == 0 else "grid"
            st.session_state["sf_result"] = compute_properties(
                tuple(X_values), tuple(T_values), tuple(m_values) if m_values else None,
                material, tuple(selected_props), mode, rhoT=rhoT, branch=branch)
            st.session_state.update(sf_ndim=ndim, sf_material=material, sf_props=selected_props,
                                    sf_range_vars=range_vars, sf_const_vars=const_vars,
                                    sf_logs=logs, sf_rhoT=rhoT)
        except Exception as e:
            st.error(f"Computation error: {e}")
            st.session_state.pop("sf_result", None)
    elif compute_btn and not selected_props:
        st.warning("Select at least one property to compute.")
    elif compute_btn and ndim > 2:
        st.error("Cannot compute with 3 independent ranges. Set one variable to a single value.")
    elif compute_btn and too_big:
        st.error(f"Too many points for (ρ, T) input with a Gibbs phase ({npoints} > "
                 f"{_GIBBS_RHOT_MAX}).")

    if "sf_result" in st.session_state:
        _show_results()


# ─────────────────────────────────────────────────────────────────────────────
# Results
# ─────────────────────────────────────────────────────────────────────────────
def _show_results():
    ss = st.session_state
    result, res_ndim, res_material = ss["sf_result"], ss["sf_ndim"], ss["sf_material"]
    res_props, res_range_vars, res_const_vars = ss["sf_props"], ss["sf_range_vars"], ss["sf_const_vars"]
    logs, res_rhoT = ss.get("sf_logs", {}), ss.get("sf_rhoT", False)
    res_avail = available_properties(res_material, res_rhoT)

    if res_ndim == 0:
        rows = []
        for prop in res_props:
            if prop in result:
                name, unit = res_avail[prop]
                rows.append({"Property": prop, "Description": name,
                             "Value": f"{result[prop].flat[0]:.6g}", "Unit": unit})
        df = pd.DataFrame(rows)
        st.subheader(f"Results — {label(res_material)}")
        st.dataframe(df, use_container_width=True, hide_index=True)
        if is_helmholtz(res_material) and "rho" in result and not np.isfinite(result["rho"]).all():
            st.caption("NaN: no fluid state of the chosen branch at this (P, T), or outside the "
                       "surface.")
        st.download_button("Download CSV", make_csv(df), file_name="seafreeze_single_point.csv",
                           mime="text/csv")
    elif res_ndim == 1:
        _show_1d(result, res_material, res_props, res_range_vars, res_const_vars, logs, res_avail)
    elif res_ndim == 2:
        _show_2d(result, res_material, res_props, res_range_vars, res_const_vars, logs, res_avail,
                 res_rhoT)


def _show_1d(result, material, props, range_vars, const_vars, logs, avail):
    var_name, var_vals, var_unit = range_vars[0]
    st.subheader(f"Results — {label(material)}")
    log_y = st.checkbox("Log y-axis", value=False, key="log_y_1d",
                        help="Log scale for the property axis (positive values only).")
    tabs = st.tabs(props) if len(props) > 1 else [st.container()]
    export = {_axis_title(var_name, var_unit): var_vals}
    for i, prop in enumerate(props):
        if prop not in result:
            continue
        name, unit = avail[prop]
        vals = np.squeeze(result[prop]).flatten()
        if len(vals) != len(var_vals):
            vals = vals[:len(var_vals)]
        export[f"{prop} ({unit})"] = vals
        with tabs[i]:
            fig = go.Figure(go.Scatter(x=var_vals, y=vals, mode="lines", name=prop,
                                       line=dict(width=2)))
            fig.update_layout(
                xaxis=dict(title=_axis_title(var_name, var_unit),
                           type="log" if logs.get(var_name) else "linear", exponentformat="power"),
                yaxis=dict(title=f"{prop} ({unit})", type="log" if log_y else "linear",
                           exponentformat="power"),
                title=title_with_fixed(name, const_vars), template="plotly_white", height=500)
            st.plotly_chart(fig, use_container_width=True)
    st.download_button("Download CSV", make_csv(pd.DataFrame(export)),
                       file_name=f"seafreeze_{material}_1D.csv", mime="text/csv")


def _fluid_overlay(material, rhoT):
    """Full-diagram overlay for a Helmholtz fluid map (precomputed default diagram)."""
    from core import diagram_plot as DP
    from core import diagrams as D
    coords = "rhoT" if rhoT else "PT"
    d = D.get_diagram(coords, D.window_items(D.DEFAULT_WINDOW[coords]), D.DEFAULT_RESOLUTION)
    return DP.overlay_traces(d, ("boundaries", "saturation", "critical", "triple"), halo=True)


def _show_2d(result, material, props, range_vars, const_vars, logs, avail, rhoT):
    var1_name, var1_vals, var1_unit = range_vars[0]
    var2_name, var2_vals, var2_unit = range_vars[1]
    st.subheader(f"Results — {label(material)}")

    valid_props = [p for p in props if p in result and np.squeeze(result[p]).ndim == 2]
    flat_options = [s for _, items in categorized_properties(material, rhoT)
                    for s, _, _ in items if s in valid_props]
    prop_fmt = {s: f"{s} — {avail[s][0]} ({avail[s][1]})" for s in flat_options}

    c1, c2, c3 = st.columns([2, 2, 1])
    with c1:
        prop = st.selectbox("Property", flat_options, format_func=lambda k: prop_fmt[k],
                            key="prop_2d")
    with c2:
        view_mode = st.radio("View", ["Heatmap", "Heatmap + Isocontours", "3-D Surface"],
                             horizontal=True, key="view_mode")
    with c3:
        colorscale = st.selectbox(
            "Color scale", ["Viridis", "Plasma", "Inferno", "Magma", "Cividis", "Turbo", "Hot",
                            "YlOrRd", "YlGnBu", "RdBu", "RdYlBu", "Spectral", "Jet", "Rainbow",
                            "Portland"], key="colorscale")
    n_contours = None
    if view_mode == "Heatmap + Isocontours":
        n_contours = st.slider("Number of contour lines", 5, 30, 15, key="n_contours")

    # stability boundaries: Gibbs phases in (P, T); water_Brown2026 in (P, T) or (rho, T)
    axis_names = {var1_name, var2_name}
    helm = is_helmholtz(material)
    can_show = axis_names == {"P", "T"} or (helm and axis_names == {"rho", "T"})
    show_boundaries, phase_bounds = False, []
    if can_show:
        show_boundaries = st.checkbox(
            "Show stability field boundaries" if not helm
            else "Show the phase diagram (boundaries, saturation, critical and triple points)",
            key="show_boundaries")
        if show_boundaries and not helm:
            phase_bounds = get_stability_boundaries(material)

    xtype = "log" if logs.get(var1_name) else "linear"
    ytype = "log" if logs.get(var2_name) else "linear"
    if prop and prop in result:
        name, unit = avail[prop]
        data_2d = np.squeeze(result[prop])
        hover = (f"{_VAR_LABEL[var1_name]}: %{{x:.4g}} {var1_unit}<br>"
                 f"{_VAR_LABEL[var2_name]}: %{{y:.4g}} {var2_unit}<br>"
                 f"{prop}: %{{z:.5g}} {unit}<extra></extra>")
        layout2d = dict(
            xaxis=dict(title=_axis_title(var1_name, var1_unit), type=xtype, exponentformat="power",
                       range=_axis_range(var1_vals, xtype)),
            yaxis=dict(title=_axis_title(var2_name, var2_unit), type=ytype, exponentformat="power",
                       range=_axis_range(var2_vals, ytype)),
            title=title_with_fixed(name, const_vars), template="plotly_white", height=600,
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0))

        if view_mode in ("Heatmap", "Heatmap + Isocontours"):
            fig = go.Figure()
            fig.add_trace(go.Heatmap(z=data_2d.T, x=var1_vals, y=var2_vals, colorscale=colorscale,
                                     colorbar=dict(title=f"{prop}<br>({unit})"),
                                     hovertemplate=hover))
            if view_mode == "Heatmap + Isocontours":
                fig.add_trace(go.Contour(
                    z=data_2d.T, x=var1_vals, y=var2_vals,
                    contours=dict(coloring="none", showlabels=True,
                                  labelfont=dict(size=11, color="white")),
                    ncontours=n_contours, showscale=False, line=dict(width=1.5, color="white"),
                    hoverinfo="skip"))
            if show_boundaries:
                if helm:
                    for tr in _fluid_overlay(material, rhoT):
                        fig.add_trace(tr)
                else:
                    add_phase_boundaries(fig, phase_bounds, var1_name, var2_name)
            fig.update_layout(**layout2d)
            st.plotly_chart(fig, use_container_width=True)
        else:  # 3-D Surface
            fig = go.Figure(go.Surface(z=data_2d.T, x=var1_vals, y=var2_vals, colorscale=colorscale,
                                       colorbar=dict(title=f"{prop}<br>({unit})")))
            if show_boundaries:
                _boundaries_3d(fig, material, rhoT, phase_bounds, var1_name, var1_vals, var2_vals,
                               data_2d)
            fig.update_layout(
                scene=dict(xaxis=dict(title=_axis_title(var1_name, var1_unit), type=xtype),
                           yaxis=dict(title=_axis_title(var2_name, var2_unit), type=ytype),
                           zaxis_title=f"{prop} ({unit})", domain=dict(x=[0, 0.85])),
                title=title_with_fixed(name, const_vars), template="plotly_white", height=700,
                legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0))
            st.plotly_chart(fig, use_container_width=True)

    # Export
    v1, v2 = np.meshgrid(var1_vals, var2_vals, indexing="ij")
    export = {_axis_title(var1_name, var1_unit): v1.flatten(),
              _axis_title(var2_name, var2_unit): v2.flatten()}
    for p in props:
        if p in result and np.size(result[p]) == v1.size:
            export[f"{p} ({avail[p][1]})"] = np.asarray(result[p]).flatten()
    st.download_button("Download CSV", make_csv(pd.DataFrame(export)),
                       file_name=f"seafreeze_{material}_2D.csv", mime="text/csv")


def _axis_range(vals, axtype):
    lo, hi = float(np.min(vals)), float(np.max(vals))
    return [np.log10(lo), np.log10(hi)] if axtype == "log" else [lo, hi]


def _boundaries_3d(fig, material, rhoT, phase_bounds, var1_name, var1_vals, var2_vals, data_2d):
    """Boundary lines on the surface (Gibbs phases) or on its floor (water_Brown2026)."""
    if is_helmholtz(material):
        zf = float(np.nanmin(data_2d)) if np.isfinite(data_2d).any() else 0.0
        x0, x1 = min(var1_vals), max(var1_vals)
        y0, y1 = min(var2_vals), max(var2_vals)
        first = True
        for tr in _fluid_overlay(material, rhoT):
            if tr.x is None or tr.line is None or tr.line.color == "white" or tr.mode != "lines":
                continue
            x, y = np.asarray(tr.x, float), np.asarray(tr.y, float)
            keep = (x >= x0) & (x <= x1) & (y >= y0) & (y <= y1)
            if keep.sum() < 2:
                continue
            fig.add_trace(go.Scatter3d(x=np.where(keep, x, np.nan), y=np.where(keep, y, np.nan),
                                       z=np.full(x.size, zf), mode="lines",
                                       line=dict(color="#0b0b0b", width=4),
                                       name="phase diagram (floor)", showlegend=first,
                                       legendgroup="floor", hoverinfo="skip"))
            first = False
        return
    from scipy.interpolate import RegularGridInterpolator
    interp = RegularGridInterpolator((np.array(var1_vals), np.array(var2_vals)), data_2d,
                                     bounds_error=False, fill_value=None)
    for j, (mat, other, bP, bT) in enumerate(phase_bounds):
        bx, by = (bP, bT) if var1_name == "P" else (bT, bP)
        bz = interp(np.column_stack([bx, by]))
        fig.add_trace(go.Scatter3d(x=bx, y=by, z=bz, mode="lines",
                                   line=dict(color=BOUNDARY_COLORS[j % len(BOUNDARY_COLORS)], width=5),
                                   name=f"{mat}–{other}"))
