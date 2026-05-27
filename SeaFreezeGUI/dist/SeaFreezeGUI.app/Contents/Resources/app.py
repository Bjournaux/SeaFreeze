"""SeaFreeze GUI — Streamlit application.

Property Calculator + Phase Diagram viewer.
"""

import io
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

import sys, os
# Ensure core/ is importable — works from source and PyInstaller bundle
_app_dir = os.path.dirname(os.path.abspath(__file__))
if _app_dir not in sys.path:
    sys.path.insert(0, _app_dir)
# PyInstaller bundles files under sys._MEIPASS
if hasattr(sys, "_MEIPASS") and sys._MEIPASS not in sys.path:
    sys.path.insert(0, sys._MEIPASS)

from core.constants import (
    ALL_MATERIALS, MATERIAL_LABELS, MATERIAL_SHORT_LABELS,
    available_properties, categorized_properties, is_nacl, is_solid,
)
from core.compute import (
    compute_properties, get_phase_range, get_stability_boundaries,
    get_phase_line,
)

# ─────────────────────────────────────────────────────────────────────────────
# Page config
# ─────────────────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="SeaFreeze",
    page_icon="❄️",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ─────────────────────────────────────────────────────────────────────────────
# Asset path resolution (works for both `streamlit run` and PyInstaller bundle)
# ─────────────────────────────────────────────────────────────────────────────
_ASSETS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets")

def _asset(name):
    if hasattr(sys, "_MEIPASS"):
        return os.path.join(sys._MEIPASS, "assets", name)
    return os.path.join(_ASSETS_DIR, name)


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _label(material):
    """Full citation label — used in selectors and page headers."""
    return MATERIAL_LABELS.get(material, material)

def _short_label(material):
    """Short label — used in plot legends to keep trace names compact."""
    return MATERIAL_SHORT_LABELS.get(material, material)

def _fmt_range(lo, hi, unit):
    return f"{lo:.1f} – {hi:.1f} {unit}"

def _make_csv(df: pd.DataFrame) -> str:
    return df.to_csv(index=False)

def _title_with_fixed(name, const_vars):
    """Build plot title appending fixed variable values."""
    if not const_vars:
        return name
    fixed = ", ".join(f"{v[0]} = {v[1]:.4g} {v[2]}" for v in const_vars)
    return f"{name}<br><sup>{fixed}</sup>"


_BOUNDARY_COLORS = [
    "#FF6347", "#FFD700", "#00FFFF", "#FF69B4",
    "#7FFF00", "#FF8C00", "#DA70D6", "#1E90FF",
]


def _add_phase_boundaries(fig, phase_bounds, var1_name, var2_name, mode="2d"):
    """Overlay phase boundary lines on a Plotly figure."""
    for i, (mat, other, bP, bT) in enumerate(phase_bounds):
        if var1_name == "P":
            x_vals, y_vals = bP, bT
        else:
            x_vals, y_vals = bT, bP
        color = _BOUNDARY_COLORS[i % len(_BOUNDARY_COLORS)]
        label = f"{mat}–{other}"
        if mode == "3d":
            fig.add_trace(go.Scatter3d(
                x=x_vals, y=y_vals,
                z=[None] * len(x_vals),
                mode="lines",
                line=dict(color=color, width=4),
                name=label,
                showlegend=True,
            ))
        else:
            fig.add_trace(go.Scatter(
                x=x_vals, y=y_vals,
                mode="lines",
                line=dict(color=color, width=2.5, dash="dot"),
                name=label,
                showlegend=True,
            ))


# ─────────────────────────────────────────────────────────────────────────────
# Sidebar branding — shown on every page
# ─────────────────────────────────────────────────────────────────────────────
with st.sidebar:
    logo_path = _asset("logo.png")
    if os.path.exists(logo_path):
        st.image(logo_path, width=230)
    st.caption(
        "Developed by **Baptiste Journaux**  \n"
        "University of Washington"
    )
    st.divider()


# ─────────────────────────────────────────────────────────────────────────────
# Page selector (top of main area)
# ─────────────────────────────────────────────────────────────────────────────
page = st.pills("Navigation",
                ["Property Calculator", "Phase Diagram", "About"],
                default="Property Calculator",
                key="page_select", label_visibility="collapsed")


# ═════════════════════════════════════════════════════════════════════════════
# PAGE 1 — Property Calculator
# ═════════════════════════════════════════════════════════════════════════════
def page_property_calculator():
    st.title("SeaFreeze — Property Calculator")

    with st.sidebar:
        st.header("Material")

        material = st.selectbox(
            "Phase / Material",
            ALL_MATERIALS,
            format_func=_label,
            key="pc_material",
        )

        nacl = is_nacl(material)
        solid = is_solid(material)

        P_range, T_range, m_range = get_phase_range(material)

        st.divider()
        st.header("Coordinates")

        st.caption(f"**P**: {_fmt_range(*P_range, 'MPa')}  \n"
                   f"**T**: {_fmt_range(*T_range, 'K')}"
                   + (f"  \n**m**: {_fmt_range(*m_range, 'mol/kg')}" if nacl else ""))

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
            st.subheader(f"Results — {_label(res_material)}")
            st.dataframe(df, use_container_width=True, hide_index=True)
            st.download_button("Download CSV", _make_csv(df),
                               file_name="seafreeze_single_point.csv",
                               mime="text/csv")

        elif res_ndim == 1:
            var_name, var_vals, var_unit = res_range_vars[0]
            st.subheader(f"Results — {_label(res_material)}")

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
                        title=_title_with_fixed(name, res_const_vars),
                        template="plotly_white",
                        height=500,
                    )
                    st.plotly_chart(fig, use_container_width=True)

            export_df = pd.DataFrame(export_data)
            st.download_button("Download CSV", _make_csv(export_df),
                               file_name=f"seafreeze_{res_material}_1D.csv",
                               mime="text/csv")

        elif res_ndim == 2:
            var1_name, var1_vals, var1_unit = res_range_vars[0]
            var2_name, var2_vals, var2_unit = res_range_vars[1]

            st.subheader(f"Results — {_label(res_material)}")

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
                        _add_phase_boundaries(fig, phase_bounds, var1_name, var2_name)
                    fig.update_layout(
                        xaxis_title=f"{var1_name} ({var1_unit})",
                        yaxis_title=f"{var2_name} ({var2_unit})",
                        title=_title_with_fixed(name, res_const_vars),
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
                        _add_phase_boundaries(fig, phase_bounds, var1_name, var2_name)
                    fig.update_layout(
                        xaxis_title=f"{var1_name} ({var1_unit})",
                        yaxis_title=f"{var2_name} ({var2_unit})",
                        title=_title_with_fixed(name, res_const_vars),
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
                            color = _BOUNDARY_COLORS[j % len(_BOUNDARY_COLORS)]
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
                        title=_title_with_fixed(name, res_const_vars),
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
            st.download_button("Download CSV", _make_csv(export_df),
                               file_name=f"seafreeze_{res_material}_2D.csv",
                               mime="text/csv")


# ═════════════════════════════════════════════════════════════════════════════
# PAGE 2 — Phase Diagram
# ═════════════════════════════════════════════════════════════════════════════

# Colors for each phase — consistent across the page
_PHASE_COLORS = {
    "Ih":           "#1f77b4",
    "II":           "#ff7f0e",
    "III":          "#2ca02c",
    "V":            "#d62728",
    "VI":           "#9467bd",
    "VII_X_French": "#8c564b",
    "water1":       "#17becf",
}

_NACL_COLORS = [
    "#e377c2", "#bcbd22", "#7f7f7f", "#aec7e8",
    "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
]

# All pure-water phase pairs (no NaCl)
_PURE_PAIRS = [
    ("Ih", "water1"), ("Ih", "II"), ("Ih", "III"),
    ("II", "III"), ("II", "V"), ("II", "VI"),
    ("III", "V"), ("III", "water1"),
    ("V", "water1"), ("V", "VI"),
    ("VI", "water1"),
]

# NaCl melting pairs
_NACL_PAIRS = [
    ("Ih", "NaClaq"), ("II", "NaClaq"), ("III", "NaClaq"),
    ("V", "NaClaq"), ("VI", "NaClaq"),
]

# All phases that can be toggled
_ALL_PHASES = ["Ih", "II", "III", "V", "VI", "water1"]


def page_phase_diagram():
    st.title("SeaFreeze — Phase Diagram")

    with st.sidebar:
        st.header("Phase boundaries")
        st.caption("Select phases to show their stability boundaries. "
                   "All boundaries between checked phases are drawn.")

        checked_phases = []
        cols = st.columns(3)
        for i, phase in enumerate(_ALL_PHASES):
            with cols[i % 3]:
                if st.checkbox(_short_label(phase), value=(phase in ("Ih", "water1")),
                               key=f"pd_{phase}"):
                    checked_phases.append(phase)

        st.divider()
        st.header("Segment")
        segment = st.radio("Segment", ["stable", "all", "metastable"],
                            horizontal=True, key="pd_segment",
                            label_visibility="collapsed")

        st.divider()
        st.header("NaCl(aq) melting curves")
        show_nacl = st.checkbox("Show NaCl(aq) curves", key="pd_show_nacl")
        nacl_molalities = []
        if show_nacl:
            nacl_m_input = st.text_input(
                "Molalities (comma-separated)",
                value="1, 2, 3",
                key="pd_nacl_m",
            )
            try:
                nacl_molalities = [float(x.strip()) for x in nacl_m_input.split(",")
                                   if x.strip()]
            except ValueError:
                st.error("Enter molalities as comma-separated numbers.")
                nacl_molalities = []

        st.divider()
        st.header("Axis ranges")
        P_lo = st.number_input("P min (MPa)", value=0.0, step=50.0, key="pd_Plo")
        P_hi = st.number_input("P max (MPa)", value=2500.0, step=50.0, key="pd_Phi")
        T_lo = st.number_input("T min (K)", value=200.0, step=10.0, key="pd_Tlo")
        T_hi = st.number_input("T max (K)", value=400.0, step=10.0, key="pd_Thi")

    # ── Build plot ────────────────────────────────────────────────────────
    fig = go.Figure()
    traces_added = 0

    # Pure-water phase boundaries between checked phases
    for matA, matB in _PURE_PAIRS:
        if matA not in checked_phases or matB not in checked_phases:
            continue
        try:
            res = get_phase_line(matA, matB, segment=segment)
            if res is None:
                continue
            P_line, T_line = res
            if len(P_line) == 0:
                continue

            dash = "solid" if segment != "metastable" else "dash"
            color = _PHASE_COLORS.get(matA, "#888888")
            fig.add_trace(go.Scatter(
                x=P_line, y=T_line,
                mode="lines",
                line=dict(color=color, width=2.5, dash=dash),
                name=f"{_short_label(matA)} – {_short_label(matB)}",
                hovertemplate=(
                    "P: %{x:.1f} MPa<br>"
                    "T: %{y:.1f} K"
                    "<extra></extra>"
                ),
            ))
            traces_added += 1
        except Exception:
            pass

    # NaCl melting curves
    if show_nacl and nacl_molalities:
        for matA, _ in _NACL_PAIRS:
            if matA not in checked_phases:
                continue
            for j, m_val in enumerate(nacl_molalities):
                try:
                    res = get_phase_line(matA, "NaClaq", segment=segment, m=m_val)
                    if res is None:
                        continue
                    P_line, T_line = res
                    if len(P_line) == 0:
                        continue
                    color = _NACL_COLORS[j % len(_NACL_COLORS)]
                    fig.add_trace(go.Scatter(
                        x=P_line, y=T_line,
                        mode="lines",
                        line=dict(color=color, width=2, dash="dot"),
                        name=f"{_short_label(matA)} – NaCl(aq) m={m_val}",
                        hovertemplate=(
                            f"m = {m_val} mol/kg<br>"
                            "T: %{x:.1f} K<br>"
                            "P: %{y:.1f} MPa"
                            "<extra></extra>"
                        ),
                    ))
                    traces_added += 1
                except Exception:
                    pass

    if traces_added == 0:
        st.info("Select at least two phases to display phase boundaries.")
    else:
        fig.update_layout(
            xaxis_title="Pressure (MPa)",
            yaxis_title="Temperature (K)",
            xaxis=dict(range=[P_lo, P_hi]),
            yaxis=dict(range=[T_lo, T_hi]),
            template="plotly_white",
            height=700,
            legend=dict(orientation="h", yanchor="bottom",
                        y=1.02, xanchor="left", x=0),
        )
        st.plotly_chart(fig, use_container_width=True)

        # Export button
        export_rows = []
        for trace in fig.data:
            for p, t in zip(trace.x, trace.y):
                export_rows.append({
                    "Boundary": trace.name,
                    "P (MPa)": p,
                    "T (K)": t,
                })
        if export_rows:
            export_df = pd.DataFrame(export_rows)
            st.download_button("Download CSV", _make_csv(export_df),
                               file_name="seafreeze_phase_diagram.csv",
                               mime="text/csv")


# ═════════════════════════════════════════════════════════════════════════════
# PAGE 3 — About & References
# ═════════════════════════════════════════════════════════════════════════════

def page_about():
    st.title("SeaFreeze — About & References")

    # ── About ─────────────────────────────────────────────────────────────────
    st.header("About")
    st.markdown("""
    **SeaFreeze GUI** is an interactive web application for computing thermodynamic
    and elastic properties of water, ice polymorphs (Ih, II, III, V, VI, VII/X), and
    NaCl aqueous solutions, built on the
    [SeaFreeze](https://github.com/Bjournaux/SeaFreeze) library.

    SeaFreeze is based on the evaluation of Gibbs Local Basis Function (LBF)
    parametrizations for each phase, constructed to reproduce thermodynamic measurements
    across a wide range of pressures and temperatures relevant to planetary interiors and
    high-pressure geophysics.

    **Try it online:** [seafreeze.streamlit.app](https://seafreeze.streamlit.app/)
    ⚠️ *Beta — under active development.*
    """)

    # ── Supported phases & validity ranges ───────────────────────────────────
    st.header("Supported Phases & Validity Ranges")
    rows = []
    for mat in ALL_MATERIALS:
        try:
            P_rng, T_rng, m_rng = get_phase_range(mat)
            row = {
                "Material": MATERIAL_LABELS[mat],
                "P range (MPa)": f"{P_rng[0]:.0f} – {P_rng[1]:.0f}",
                "T range (K)":   f"{T_rng[0]:.0f} – {T_rng[1]:.0f}",
            }
            if m_rng is not None:
                row["m range (mol/kg)"] = f"{m_rng[0]:.1f} – {m_rng[1]:.1f}"
            rows.append(row)
        except Exception:
            pass
    st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)

    st.caption(
        "NaN values are returned outside the parametrization boundaries. "
        "Stability predictions are considered valid down to 130 K (ice VI–XV transition). "
        "The ice Ih–II transition is potentially valid down to 73.4 K."
    )

    # ── Phase diagram validation ──────────────────────────────────────────────
    st.header("Phase Diagram — Validation Against Experimental Data")
    try:
        import seafreeze as _sf_pkg
        fig_path = os.path.join(os.path.dirname(_sf_pkg.__file__),
                                "figures", "WPD_python.png")
        if os.path.exists(fig_path):
            st.image(fig_path,
                     caption="SeaFreeze computed phase boundaries vs. experimental data.",
                     use_container_width=True)
    except Exception:
        st.info("Phase diagram figure not available in this environment.")

    st.markdown("""
    The ice Ih–VII/X melting curve above the VI–VII–water triple point (~2216 MPa, 354 K)
    uses the Gibbs energy representation of French & Redmer (2015).

    For phase equilibrium calculations, **water1** (Bollengier et al., 2019) is recommended
    over water2 or IAPWS95, as the ice Gibbs parametrizations are optimized against it.
    """)

    # ── Contributors ─────────────────────────────────────────────────────────
    st.header("Contributors")
    st.markdown("""
    - **Baptiste Journaux** *(Lead)* — University of Washington, Earth & Space Sciences
    - **J. Michael Brown** — University of Washington, Earth & Space Sciences
    - **Penny Espinoza** — University of Washington, Earth & Space Sciences
    - **Ula Jones** — University of Washington, Earth & Space Sciences
    - **Erica Clinton** — University of Washington, Earth & Space Sciences
    - **Tyler Gordon** — University of Washington, Department of Astronomy
    - **Matthew J. Powell-Palm** — Texas A&M University, Mechanical Engineering
    - **Steven D. Vance** — NASA Jet Propulsion Laboratory, Caltech
    """)

    # ── References ───────────────────────────────────────────────────────────
    st.header("References")
    st.markdown("""
    1. **Journaux et al. (2020)** — Holistic Approach for Studying Planetary Hydrospheres:
       Gibbs Representation of Ices Thermodynamics, Elasticity, and the Water Phase Diagram
       to 2300 MPa. *JGR Planets* 125(1), e2019JE006176.
       [DOI: 10.1029/2019JE006176](https://doi.org/10.1029/2019JE006176)

    2. **Bollengier, Brown & Shaw (2019)** — Thermodynamics of pure liquid water to 2300 MPa
       and 500 K from a new Gibbs parametrization.
       *J. Chem. Phys.* 151, 054501.
       [DOI: 10.1063/1.5097179](https://doi.org/10.1063/1.5097179)

    3. **Brown (2018)** — Seismic wave anisotropy in the inner core and its relation to
       high-pressure solid phase transitions.
       *Fluid Phase Equilibria* 463, 18–31.
       [DOI: 10.1016/j.fluid.2018.02.001](https://doi.org/10.1016/j.fluid.2018.02.001)

    4. **Feistel & Wagner (2006)** — A New Equation of State for H₂O Ice Ih.
       *J. Phys. Chem. Ref. Data* 35, 1021–1047.

    5. **Wagner & Pruss (2002)** — The IAPWS Formulation 1995 for the Thermodynamic
       Properties of Ordinary Water Substance for General and Scientific Use.
       *J. Phys. Chem. Ref. Data* 31, 387–535.

    6. **French & Redmer (2015)** — Electronic band structure and optical properties of
       ice VII and X.
       *Physical Review B* 91, 014308.
       [DOI: 10.1103/PhysRevB.91.014308](https://doi.org/10.1103/PhysRevB.91.014308)

    7. **Brown et al. (under review)** — NaCl aqueous solution equation of state.
    """)


# ─────────────────────────────────────────────────────────────────────────────
# Route to the selected page
# ─────────────────────────────────────────────────────────────────────────────
if page == "Property Calculator":
    page_property_calculator()
elif page == "Phase Diagram":
    page_phase_diagram()
else:
    page_about()
