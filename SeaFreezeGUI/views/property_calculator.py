"""Property Calculator page — single points, 1-D sweeps and 2-D grids."""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from core.constants import (
    ALL_MATERIALS, available_properties, categorized_properties, is_nacl,
    is_solid,
)
from core.compute import (
    compute_properties, get_phase_range, get_stability_boundaries,
)
from core.ui import (
    BOUNDARY_COLORS, add_phase_boundaries, fmt_range, label, make_csv,
    title_with_fixed,
)


def render():
    st.title("SeaFreeze — Property Calculator")

    with st.sidebar:
        st.header("Material")

        material = st.selectbox(
            "Phase / Material",
            ALL_MATERIALS,
            format_func=label,
            key="pc_material",
        )

        nacl = is_nacl(material)
        solid = is_solid(material)

        P_range, T_range, m_range = get_phase_range(material)

        st.divider()
        st.header("Coordinates")

        st.caption(f"**P**: {fmt_range(*P_range, 'MPa')}  \n"
                   f"**T**: {fmt_range(*T_range, 'K')}"
                   + (f"  \n**m**: {fmt_range(*m_range, 'mol/kg')}" if nacl else ""))

        # ── Pressure ──────────────────────────────────────────────────────
        st.subheader("Pressure (MPa)")
        P_mode = st.radio("P input", ["Single value", "Range"],
                           key="P_mode", horizontal=True, label_visibility="collapsed")
        if P_mode == "Single value":
            P_val = st.number_input("P (MPa)",
                                     value=float(round((P_range[0] + P_range[1]) / 2, 1)),
                                     min_value=float(P_range[0]),
                                     max_value=float(P_range[1]),
                                     step=10.0, key="P_single")
            P_values = [P_val]
            P_is_range = False
        else:
            pc1, pc2 = st.columns(2)
            P_min = pc1.number_input("P min", value=float(P_range[0]),
                                      min_value=float(P_range[0]),
                                      max_value=float(P_range[1]),
                                      step=10.0, key="P_min")
            P_max = pc2.number_input("P max", value=float(P_range[1]),
                                      min_value=float(P_range[0]),
                                      max_value=float(P_range[1]),
                                      step=10.0, key="P_max")
            P_npts = st.number_input("P points", value=50, min_value=2,
                                      max_value=500, key="P_npts")
            P_values = np.linspace(P_min, P_max, int(P_npts)).tolist()
            P_is_range = True

        # ── Temperature ───────────────────────────────────────────────────
        st.subheader("Temperature (K)")
        T_mode = st.radio("T input", ["Single value", "Range"],
                           key="T_mode", horizontal=True, label_visibility="collapsed")
        if T_mode == "Single value":
            T_val = st.number_input("T (K)",
                                     value=float(round((T_range[0] + T_range[1]) / 2, 1)),
                                     min_value=float(T_range[0]),
                                     max_value=float(T_range[1]),
                                     step=5.0, key="T_single")
            T_values = [T_val]
            T_is_range = False
        else:
            tc1, tc2 = st.columns(2)
            T_min = tc1.number_input("T min", value=float(T_range[0]),
                                      min_value=float(T_range[0]),
                                      max_value=float(T_range[1]),
                                      step=5.0, key="T_min")
            T_max = tc2.number_input("T max", value=float(T_range[1]),
                                      min_value=float(T_range[0]),
                                      max_value=float(T_range[1]),
                                      step=5.0, key="T_max")
            T_npts = st.number_input("T points", value=50, min_value=2,
                                      max_value=500, key="T_npts")
            T_values = np.linspace(T_min, T_max, int(T_npts)).tolist()
            T_is_range = True

        # ── Molality (NaClaq only) ────────────────────────────────────────
        m_values = None
        m_is_range = False
        if nacl:
            st.subheader("Molality (mol/kg)")
            m_mode = st.radio("m input", ["Single value", "Range"],
                               key="m_mode", horizontal=True,
                               label_visibility="collapsed")
            if m_mode == "Single value":
                m_val = st.number_input("m (mol/kg)", value=1.0,
                                         min_value=float(m_range[0]),
                                         max_value=float(m_range[1]),
                                         step=0.5, key="m_single")
                m_values = [m_val]
                m_is_range = False
            else:
                mc1, mc2 = st.columns(2)
                m_min = mc1.number_input("m min", value=float(m_range[0]),
                                          min_value=float(m_range[0]),
                                          max_value=float(m_range[1]),
                                          step=0.5, key="m_min")
                m_max = mc2.number_input("m max", value=float(m_range[1]),
                                          min_value=float(m_range[0]),
                                          max_value=float(m_range[1]),
                                          step=0.5, key="m_max")
                m_npts = st.number_input("m points", value=20, min_value=2,
                                          max_value=200, key="m_npts")
                m_values = np.linspace(m_min, m_max, int(m_npts)).tolist()
                m_is_range = True

        # ── Property selection ────────────────────────────────────────────
        st.divider()
        st.header("Properties")
        avail = available_properties(material)
        prop_labels = {k: f"{k}  ({v[0]}, {v[1]})" for k, v in avail.items()}

        select_all = st.checkbox("Compute all properties", value=True,
                                  key="select_all")
        if select_all:
            selected_props = list(avail.keys())
        else:
            selected_props = st.multiselect(
                "Select properties",
                options=list(avail.keys()),
                default=["rho", "Cp", "G"],
                format_func=lambda k: prop_labels[k],
                key="prop_select",
            )

        # ── Compute button ────────────────────────────────────────────────
        st.divider()
        compute_btn = st.button("Compute", type="primary",
                                 use_container_width=True)

    # ── Determine dimensionality ──────────────────────────────────────────
    range_vars = []
    const_vars = []

    if P_is_range:
        range_vars.append(("P", P_values, "MPa"))
    else:
        const_vars.append(("P", P_values[0], "MPa"))

    if T_is_range:
        range_vars.append(("T", T_values, "K"))
    else:
        const_vars.append(("T", T_values[0], "K"))

    if nacl:
        if m_is_range:
            range_vars.append(("m", m_values, "mol/kg"))
        else:
            const_vars.append(("m", m_values[0], "mol/kg"))

    ndim = len(range_vars)

    mode_labels = {
        0: "Single point",
        1: "1-D sweep (line plots)",
        2: "2-D grid (heatmap / 3-D surface)",
    }
    if ndim <= 2:
        mode_str = f"**Mode**: {mode_labels[ndim]}"
        if ndim > 0:
            mode_str += f" — varying **{' & '.join(v[0] for v in range_vars)}**"
        if const_vars:
            mode_str += " | " + ", ".join(
                f"{v[0]} = {v[1]} {v[2]}" for v in const_vars)
        st.info(mode_str)
    else:
        st.warning(
            "3 independent ranges selected. Only 2 can vary at once for a "
            "surface plot. Set one variable to a single value.")

    # ── Compute & store ───────────────────────────────────────────────────
    if compute_btn and selected_props and ndim <= 2:
        try:
            props_tuple = tuple(selected_props)
            mode = "scatter" if ndim == 0 else "grid"
            st.session_state["sf_result"] = compute_properties(
                tuple(P_values), tuple(T_values),
                tuple(m_values) if m_values else None,
                material, props_tuple, mode,
            )
            st.session_state["sf_ndim"] = ndim
            st.session_state["sf_material"] = material
            st.session_state["sf_props"] = selected_props
            st.session_state["sf_range_vars"] = range_vars
            st.session_state["sf_const_vars"] = const_vars
        except Exception as e:
            st.error(f"Computation error: {e}")
            st.session_state.pop("sf_result", None)

    elif compute_btn and not selected_props:
        st.warning("Select at least one property to compute.")
    elif compute_btn and ndim > 2:
        st.error("Cannot compute with 3 independent ranges. "
                 "Set one variable to a single value.")

    # ── Display results from session state ────────────────────────────────
    if "sf_result" in st.session_state:
        result = st.session_state["sf_result"]
        res_ndim = st.session_state["sf_ndim"]
        res_material = st.session_state["sf_material"]
        res_props = st.session_state["sf_props"]
        res_range_vars = st.session_state["sf_range_vars"]
        res_const_vars = st.session_state["sf_const_vars"]
        res_avail = available_properties(res_material)

        if res_ndim == 0:
            rows = []
            for prop in res_props:
                if prop in result:
                    val = result[prop]
                    name, unit = res_avail[prop]
                    rows.append({
                        "Property": prop,
                        "Description": name,
                        "Value": f"{val.flat[0]:.6g}",
                        "Unit": unit,
                    })
            df = pd.DataFrame(rows)
            st.subheader(f"Results — {label(res_material)}")
            st.dataframe(df, use_container_width=True, hide_index=True)
            st.download_button("Download CSV", make_csv(df),
                               file_name="seafreeze_single_point.csv",
                               mime="text/csv")

        elif res_ndim == 1:
            var_name, var_vals, var_unit = res_range_vars[0]
            st.subheader(f"Results — {label(res_material)}")

            if len(res_props) > 1:
                tabs = st.tabs(res_props)
            else:
                tabs = [st.container()]

            export_data = {f"{var_name} ({var_unit})": var_vals}

            for i, prop in enumerate(res_props):
                if prop not in result:
                    continue
                arr = result[prop]
                name, unit = res_avail[prop]
                vals = np.squeeze(arr).flatten()
                if len(vals) != len(var_vals):
                    vals = vals[:len(var_vals)]
                export_data[f"{prop} ({unit})"] = vals

                with tabs[i]:
                    fig = go.Figure()
                    fig.add_trace(go.Scatter(
                        x=var_vals, y=vals,
                        mode="lines", name=prop,
                        line=dict(width=2),
                    ))
                    fig.update_layout(
                        xaxis_title=f"{var_name} ({var_unit})",
                        yaxis_title=f"{prop} ({unit})",
                        title=title_with_fixed(name, res_const_vars),
                        template="plotly_white",
                        height=500,
                    )
                    st.plotly_chart(fig, use_container_width=True)

            export_df = pd.DataFrame(export_data)
            st.download_button("Download CSV", make_csv(export_df),
                               file_name=f"seafreeze_{res_material}_1D.csv",
                               mime="text/csv")

        elif res_ndim == 2:
            var1_name, var1_vals, var1_unit = res_range_vars[0]
            var2_name, var2_vals, var2_unit = res_range_vars[1]

            st.subheader(f"Results — {label(res_material)}")

            valid_props = [p for p in res_props if p in result
                           and np.squeeze(result[p]).ndim == 2]
            categories = categorized_properties(res_material)
            grouped_options = []
            for cat_name, items in categories:
                cat_syms = [s for s, _, _ in items if s in valid_props]
                if cat_syms:
                    grouped_options.append((cat_name, cat_syms))
            flat_options = [s for _, syms in grouped_options for s in syms]
            prop_fmt = {s: f"{s} — {res_avail[s][0]} ({res_avail[s][1]})"
                        for s in flat_options}

            c1, c2, c3 = st.columns([2, 2, 1])
            with c1:
                prop = st.selectbox(
                    "Property", flat_options,
                    format_func=lambda k: prop_fmt[k],
                    key="prop_2d",
                )
            with c2:
                view_mode = st.radio(
                    "View",
                    ["Heatmap", "Heatmap + Isocontours", "3-D Surface"],
                    horizontal=True, key="view_mode",
                )
            with c3:
                colorscale = st.selectbox(
                    "Color scale",
                    ["Viridis", "Plasma", "Inferno", "Magma", "Cividis",
                     "Turbo", "Hot", "YlOrRd", "YlGnBu", "RdBu",
                     "RdYlBu", "Spectral", "Jet", "Rainbow", "Portland"],
                    key="colorscale",
                )

            n_contours = None
            if view_mode == "Heatmap + Isocontours":
                n_contours = st.slider("Number of contour lines",
                                        5, 30, 15, key="n_contours")

            axis_names = {var1_name, var2_name}
            can_show_boundaries = axis_names == {"P", "T"}
            show_boundaries = False
            phase_bounds = []
            if can_show_boundaries:
                show_boundaries = st.checkbox("Show stability field boundaries",
                                              key="show_boundaries")
                if show_boundaries:
                    phase_bounds = get_stability_boundaries(res_material)

            if prop and prop in result:
                arr = result[prop]
                name, unit = res_avail[prop]
                data_2d = np.squeeze(arr)

                if view_mode == "Heatmap":
                    fig = go.Figure(data=go.Heatmap(
                        z=data_2d.T,
                        x=var1_vals,
                        y=var2_vals,
                        colorscale=colorscale,
                        colorbar=dict(title=f"{prop}<br>({unit})"),
                        hovertemplate=(
                            f"{var1_name}: %{{x:.2f}} {var1_unit}<br>"
                            f"{var2_name}: %{{y:.2f}} {var2_unit}<br>"
                            f"{prop}: %{{z:.4g}} {unit}"
                            "<extra></extra>"
                        ),
                    ))
                    if show_boundaries:
                        add_phase_boundaries(fig, phase_bounds, var1_name, var2_name)
                    fig.update_layout(
                        xaxis_title=f"{var1_name} ({var1_unit})",
                        yaxis_title=f"{var2_name} ({var2_unit})",
                        title=title_with_fixed(name, res_const_vars),
                        template="plotly_white",
                        height=600,
                        legend=dict(orientation="h", yanchor="bottom",
                                    y=1.02, xanchor="left", x=0),
                    )
                    st.plotly_chart(fig, use_container_width=True)

                elif view_mode == "Heatmap + Isocontours":
                    fig = go.Figure()
                    fig.add_trace(go.Heatmap(
                        z=data_2d.T,
                        x=var1_vals,
                        y=var2_vals,
                        colorscale=colorscale,
                        colorbar=dict(title=f"{prop}<br>({unit})"),
                        hovertemplate=(
                            f"{var1_name}: %{{x:.2f}} {var1_unit}<br>"
                            f"{var2_name}: %{{y:.2f}} {var2_unit}<br>"
                            f"{prop}: %{{z:.4g}} {unit}"
                            "<extra></extra>"
                        ),
                    ))
                    fig.add_trace(go.Contour(
                        z=data_2d.T,
                        x=var1_vals,
                        y=var2_vals,
                        contours=dict(
                            coloring="none",
                            showlabels=True,
                            labelfont=dict(size=11, color="white"),
                        ),
                        ncontours=n_contours,
                        showscale=False,
                        line=dict(width=1.5, color="white"),
                        hoverinfo="skip",
                    ))
                    if show_boundaries:
                        add_phase_boundaries(fig, phase_bounds, var1_name, var2_name)
                    fig.update_layout(
                        xaxis_title=f"{var1_name} ({var1_unit})",
                        yaxis_title=f"{var2_name} ({var2_unit})",
                        title=title_with_fixed(name, res_const_vars),
                        template="plotly_white",
                        height=600,
                        legend=dict(orientation="h", yanchor="bottom",
                                    y=1.02, xanchor="left", x=0),
                    )
                    st.plotly_chart(fig, use_container_width=True)

                else:  # 3-D Surface
                    fig = go.Figure(data=go.Surface(
                        z=data_2d.T,
                        x=var1_vals,
                        y=var2_vals,
                        colorscale=colorscale,
                        colorbar=dict(title=f"{prop}<br>({unit})"),
                    ))
                    if show_boundaries:
                        from scipy.interpolate import RegularGridInterpolator
                        interp = RegularGridInterpolator(
                            (np.array(var1_vals), np.array(var2_vals)),
                            data_2d, bounds_error=False, fill_value=None,
                        )
                        for j, (mat, other, bP, bT) in enumerate(phase_bounds):
                            if var1_name == "P":
                                bx, by = bP, bT
                            else:
                                bx, by = bT, bP
                            pts = np.column_stack([bx, by])
                            bz = interp(pts)
                            color = BOUNDARY_COLORS[j % len(BOUNDARY_COLORS)]
                            fig.add_trace(go.Scatter3d(
                                x=bx, y=by, z=bz,
                                mode="lines",
                                line=dict(color=color, width=5),
                                name=f"{mat}–{other}",
                            ))
                    fig.update_layout(
                        scene=dict(
                            xaxis_title=f"{var1_name} ({var1_unit})",
                            yaxis_title=f"{var2_name} ({var2_unit})",
                            zaxis_title=f"{prop} ({unit})",
                            domain=dict(x=[0, 0.85]),
                        ),
                        title=title_with_fixed(name, res_const_vars),
                        template="plotly_white",
                        height=700,
                        legend=dict(orientation="h", yanchor="bottom",
                                    y=1.02, xanchor="left", x=0),
                    )
                    st.plotly_chart(fig, use_container_width=True)

            # Export
            v1_grid, v2_grid = np.meshgrid(var1_vals, var2_vals,
                                            indexing="ij")
            export_data = {
                f"{var1_name} ({var1_unit})": v1_grid.flatten(),
                f"{var2_name} ({var2_unit})": v2_grid.flatten(),
            }
            for prop in res_props:
                if prop in result:
                    _, unit = res_avail[prop]
                    export_data[f"{prop} ({unit})"] = result[prop].flatten()
            export_df = pd.DataFrame(export_data)
            st.download_button("Download CSV", make_csv(export_df),
                               file_name=f"seafreeze_{res_material}_2D.csv",
                               mime="text/csv")
