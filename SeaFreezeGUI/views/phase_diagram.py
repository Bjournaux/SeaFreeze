"""Phase Diagram page — ice/water boundaries, NaCl(aq) melting curves and
transition jumps on click."""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from core.compute import compute_properties, get_phase_line, get_phase_line_full
from core.ui import make_csv, short_label
from views import full_diagram

_MODES = ["Full diagram (water_Brown2026 + ices)", "Ice melting lines"]


def render():
    st.title("SeaFreeze — Phase Diagram")
    mode = st.radio("Mode", _MODES, horizontal=True, key="pd_mode",
                    label_visibility="collapsed")
    if mode == _MODES[0]:
        full_diagram.render()
    else:
        _render_melting_lines()


# Colors for each phase — consistent across the page
_PHASE_COLORS = {
    "Ih":           "#1f77b4",
    "II":           "#ff7f0e",
    "III":          "#2ca02c",
    "V":            "#d62728",
    "VI":           "#9467bd",
    "VII_X_French": "#8c564b",
    "water_Bollengier2019":       "#17becf",
}

# Phase-transition line colors: stable = black, metastable = dark grey
_STABLE_COLOR = "#000000"
_META_COLOR = "#777777"
_LABEL_COLOR = "#404040"   # dark grey for phase-region labels

_NACL_COLORS = [
    "#e377c2", "#bcbd22", "#7f7f7f", "#aec7e8",
    "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
]

# All pure-water phase pairs (no NaCl)
_PURE_PAIRS = [
    ("Ih", "water_Bollengier2019"), ("Ih", "II"), ("Ih", "III"),
    ("II", "III"), ("II", "V"), ("II", "VI"),
    ("III", "V"), ("III", "water_Bollengier2019"),
    ("V", "water_Bollengier2019"), ("V", "VI"),
    ("VI", "water_Bollengier2019"),
]

# NaCl melting pairs
_NACL_PAIRS = [
    ("Ih", "NaClaq"), ("II", "NaClaq"), ("III", "NaClaq"),
    ("V", "NaClaq"), ("VI", "NaClaq"),
]

# All phases that can be toggled
_ALL_PHASES = ["Ih", "II", "III", "V", "VI", "water_Bollengier2019"]

# Representative (P_MPa, T_K) interior points + short text for region labels.
_PHASE_LABEL = {
    "Ih":     (60.0,   235.0, "Ih"),
    "II":     (300.0,  218.0, "II"),
    "III":    (285.0,  249.0, "III"),
    "V":      (490.0,  248.0, "V"),
    "VI":     (1500.0, 290.0, "VI"),
    "water_Bollengier2019": (250.0,  335.0, "Liquid"),
}


_LIQUIDS = {"water_Bollengier2019", "water_Brown2026", "NaClaq"}


def _add_boundary_traces(fig, P, T, stable, color, name, custom_row,
                         segment, base_width=1.75):
    """Add the boundary line(s) for one phase pair.

    - segment == 'stable': single solid line (P/T already stable-only).
    - segment == 'all': full curve drawn as a thin dashed line (the
      metastable extension) with the stable portion overlaid as a thick
      solid line on top — so there are no gaps at the stable/meta junction.
    Returns the number of traces added.
    """
    P = np.asarray(P, dtype=float)
    T = np.asarray(T, dtype=float)
    cd = [custom_row] * len(P)
    hover = "P: %{x:.1f} MPa<br>T: %{y:.1f} K<extra></extra>"
    meta_width = max(0.6, base_width * 0.6)
    n = 0

    if segment != "all":
        fig.add_trace(go.Scatter(
            x=P, y=T, mode="lines",
            line=dict(color=color, width=base_width, dash="solid"),
            name=name, customdata=cd, hovertemplate=hover,
        ))
        return 1

    stab = (np.asarray(stable, dtype=bool) if stable is not None
            else np.ones(len(P), dtype=bool))
    has_stable = bool(np.any(stab))
    has_meta = bool(not np.all(stab))

    # Full curve as a thin dashed line (continuous underlay = no gaps).
    fig.add_trace(go.Scatter(
        x=P, y=T, mode="lines",
        line=dict(color=_META_COLOR, width=meta_width, dash="dash"),
        opacity=0.85,
        name=(name if not has_stable else f"{name} (meta)"),
        showlegend=(not has_stable),
        customdata=cd, hovertemplate=hover,
    ))
    n += 1

    # Stable portion overlaid as a solid line.
    if has_stable:
        Ps = np.where(stab, P, np.nan)
        Ts = np.where(stab, T, np.nan)
        fig.add_trace(go.Scatter(
            x=Ps, y=Ts, mode="lines",
            line=dict(color=color, width=base_width, dash="solid"),
            name=name, customdata=cd, hovertemplate=hover,
        ))
        n += 1
    return n


def _render_transition_jump(event):
    """Compute and display ΔV / ΔS / ΔH (final − initial) for the clicked point.

    Convention:
      - Melting (ice ↔ liquid): final = liquid, initial = ice.
      - Solid–solid (ice ↔ ice): final = denser (higher-P) phase.
    All values are specific (per-kg), so V = 1/rho.
    """
    try:
        pts = (event.selection or {}).get("points", []) if event else []
    except Exception:
        pts = []
    if not pts:
        return

    pt = pts[0]
    try:
        P_clk = float(pt["x"])
        T_clk = float(pt["y"])
        matA, matB, kind, m_str = pt["customdata"]
        if kind == "tp":      # triple-point marker — no transition to compute
            return
        m_val = float(m_str) if m_str else None
    except (KeyError, ValueError, TypeError):
        st.warning("Could not read the clicked point — try clicking directly "
                   "on a boundary curve.")
        return

    def _props(mat):
        m_arg = (m_val,) if mat == "NaClaq" else None
        res = compute_properties((P_clk,), (T_clk,), m_arg, mat,
                                 ("rho", "S", "H"), "scatter")
        return (float(np.ravel(res["rho"])[0]),
                float(np.ravel(res["S"])[0]),
                float(np.ravel(res["H"])[0]))

    try:
        rhoA, SA, HA = _props(matA)
        rhoB, SB, HB = _props(matB)
    except Exception as e:
        st.warning(f"Could not compute properties at this point: {e}")
        return

    if any(not np.isfinite(v) for v in (rhoA, SA, HA, rhoB, SB, HB)):
        st.warning("One of the phases is outside its valid domain at this "
                   "point — no transition jump available here.")
        return

    # Decide final vs initial
    a_liq = matA in _LIQUIDS
    b_liq = matB in _LIQUIDS
    if a_liq != b_liq:           # melting: liquid is final
        if a_liq:
            final, initial = (matA, rhoA, SA, HA), (matB, rhoB, SB, HB)
        else:
            final, initial = (matB, rhoB, SB, HB), (matA, rhoA, SA, HA)
    else:                        # solid–solid: denser phase is final
        if rhoA >= rhoB:
            final, initial = (matA, rhoA, SA, HA), (matB, rhoB, SB, HB)
        else:
            final, initial = (matB, rhoB, SB, HB), (matA, rhoA, SA, HA)

    matF, rhoF, SF_, HF = final
    matI, rhoI, SI, HI = initial
    dV = 1.0 / rhoF - 1.0 / rhoI
    dS = SF_ - SI
    dH = HF - HI

    st.divider()
    st.subheader("Transition jump (final − initial)")
    st.markdown(
        f"At **P = {P_clk:.1f} MPa**, **T = {T_clk:.2f} K** &nbsp;—&nbsp; "
        f"**{short_label(matI)} → {short_label(matF)}**")
    c1, c2, c3 = st.columns(3)
    c1.metric("ΔV (m³/kg)", f"{dV:.4e}")
    c2.metric("ΔS (J/kg/K)", f"{dS:.2f}")
    c3.metric("ΔH (J/kg)", f"{dH:.4e}", help=f"{dH / 1000.0:.2f} kJ/kg")
    if kind == "nacl":
        st.caption(
            f"NaCl(aq) end-member values: pure ice {short_label(matI)} vs "
            f"NaCl(aq) solution at m = {m_val} mol/kg. Specific (per-kg) "
            "quantities for each end-member at the clicked (P, T).")


def _render_melting_lines():
    with st.sidebar:
        st.header("Liquid")
        liquid = st.radio("Liquid", ["water_Bollengier2019", "water_Brown2026"], horizontal=True, key="pd_liquid",
                          format_func=short_label, label_visibility="collapsed",
                          help="Liquid used for the melting curves. water_Brown2026 (Helmholtz, beta) "
                               "reproduces the water_Bollengier2019 melting curves within 0.06 K to 632 MPa.")
        pure_pairs = [(a, liquid if b == "water_Bollengier2019" else b) for a, b in _PURE_PAIRS]
        all_phases = [liquid if p == "water_Bollengier2019" else p for p in _ALL_PHASES]

        st.divider()
        st.header("Phase boundaries")
        st.caption("Select phases to show their stability boundaries. "
                   "All boundaries between checked phases are drawn.")

        checked_phases = []
        cols = st.columns(3)
        for i, phase in enumerate(all_phases):
            key = "pd_liquid_on" if phase == liquid else f"pd_{phase}"
            with cols[i % 3]:
                if st.checkbox(short_label(phase), value=True, key=key):
                    checked_phases.append(phase)

        st.divider()
        st.header("Segment")
        segment = st.radio("Segment", ["stable", "all"],
                            horizontal=True, key="pd_segment",
                            format_func=lambda s: ("stable + metastable"
                                                   if s == "all" else "stable"),
                            label_visibility="collapsed")
        show_triple_points = st.checkbox("Show triple points",
                                         value=False, key="pd_show_tp")

        st.divider()
        st.header("Style")
        line_width = st.slider("Line width", min_value=0.5, max_value=4.0,
                               value=1.75, step=0.25, key="pd_lw")
        show_labels = st.checkbox("Show phase labels", value=True,
                                  key="pd_labels")

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
        st.header("Axes")
        c1, c2 = st.columns(2)
        log_P = c1.checkbox("log P", value=False, key="pd_logP")
        log_T = c2.checkbox("log T", value=False, key="pd_logT")
        st.header("Axis ranges")
        auto_axes = st.checkbox("Auto-fit axes to selected phases",
                                value=True, key="pd_autoaxes")
        P_lo = P_hi = T_lo = T_hi = None
        if not auto_axes:
            P_lo = st.number_input("P min (MPa)", value=0.0, step=50.0, key="pd_Plo")
            P_hi = st.number_input("P max (MPa)", value=2500.0, step=50.0, key="pd_Phi")
            T_lo = st.number_input("T min (K)", value=200.0, step=10.0, key="pd_Tlo")
            T_hi = st.number_input("T max (K)", value=400.0, step=10.0, key="pd_Thi")

    # ── Build plot ────────────────────────────────────────────────────────
    fig = go.Figure()
    traces_added = 0
    all_P, all_T = [], []
    triple_pts = {}   # (P_round, T_round) -> (P, T), deduped

    # Pure-water phase boundaries between checked phases
    for matA, matB in pure_pairs:
        if matA not in checked_phases or matB not in checked_phases:
            continue
        try:
            res = get_phase_line_full(matA, matB, segment=segment)
            if res is None:
                continue
            P_line, T_line, stable, tps = res
            if len(P_line) == 0:
                continue

            name = f"{short_label(matA)} – {short_label(matB)}"
            traces_added += _add_boundary_traces(
                fig, P_line, T_line, stable, _STABLE_COLOR, name,
                [matA, matB, "pure", ""], segment, base_width=line_width)
            all_P.extend(P_line); all_T.extend(T_line)

            # Collect triple points (only those within the plotted curve span)
            if tps is not None and len(tps) > 0:
                for p_tp, t_tp in np.asarray(tps, dtype=float):
                    triple_pts[(round(p_tp, 2), round(t_tp, 2))] = (p_tp, t_tp)
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
                        line=dict(color=color, width=line_width, dash="dot"),
                        name=f"{short_label(matA)} – NaCl(aq) m={m_val}",
                        customdata=[[matA, "NaClaq", "nacl", str(m_val)]] * len(P_line),
                        hovertemplate=(
                            f"m = {m_val} mol/kg<br>"
                            "P: %{x:.1f} MPa<br>"
                            "T: %{y:.1f} K"
                            "<extra></extra>"
                        ),
                    ))
                    traces_added += 1
                    all_P.extend(P_line); all_T.extend(T_line)
                except Exception:
                    pass

    # Triple points — markers shared by the drawn pure-phase boundaries
    if show_triple_points and triple_pts:
        tp_P = [v[0] for v in triple_pts.values()]
        tp_T = [v[1] for v in triple_pts.values()]
        fig.add_trace(go.Scatter(
            x=tp_P, y=tp_T, mode="markers",
            marker=dict(symbol="triangle-up", size=11, color="black",
                        line=dict(width=0.5, color="white")),
            name="Triple points",
            customdata=[["", "", "tp", ""]] * len(tp_P),
            hovertemplate=(
                "Triple point<br>P: %{x:.2f} MPa<br>T: %{y:.2f} K"
                "<extra></extra>"
            ),
        ))

    if traces_added == 0:
        st.info("Select at least two phases to display phase boundaries.")
    else:
        if auto_axes and all_P:
            def _pad(lo, hi, frac=0.05, floor=None):
                span = max(hi - lo, 1e-9)
                pad = max(span * frac, 1e-9)
                lo2 = lo - pad
                if floor is not None:
                    lo2 = max(lo2, floor)
                return lo2, hi + pad
            P_lo, P_hi = _pad(min(all_P), max(all_P), floor=0.0)
            T_lo, T_hi = _pad(min(all_T), max(all_T))
        if log_P:
            pos = [p for p in all_P if p > 0]
            P_lo = P_lo if P_lo and P_lo > 0 else (min(pos) if pos else 0.1)
        if log_T:
            T_lo = T_lo if T_lo and T_lo > 0 else 1.0
        rng = lambda lo, hi, log: [np.log10(lo), np.log10(hi)] if log else [lo, hi]
        fig.update_layout(
            xaxis_title="Pressure (MPa)",
            yaxis_title="Temperature (K)",
            xaxis=dict(range=rng(P_lo, P_hi, log_P), type="log" if log_P else "linear"),
            yaxis=dict(range=rng(T_lo, T_hi, log_T), type="log" if log_T else "linear"),
            template="plotly_white",
            height=700,
            clickmode="event+select",
            dragmode=False,
            legend=dict(orientation="h", yanchor="bottom",
                        y=1.02, xanchor="left", x=0),
        )

        # Phase region labels (only for checked phases inside the visible axes)
        if show_labels:
            for phase in checked_phases:
                key = "water_Bollengier2019" if phase == liquid else phase
                if key not in _PHASE_LABEL:
                    continue
                Px, Tx, txt = _PHASE_LABEL[key]
                if not (P_lo <= Px <= P_hi and T_lo <= Tx <= T_hi):
                    continue
                fig.add_annotation(
                    x=np.log10(Px) if log_P else Px, y=np.log10(Tx) if log_T else Tx,
                    text=f"<b>{txt}</b>", showarrow=False,
                    font=dict(size=20, color=_LABEL_COLOR),
                )

        st.caption("Click a point on any boundary curve to compute the "
                   "thermodynamic jump (ΔV, ΔS, ΔH) across the transition.")
        event = st.plotly_chart(
            fig, use_container_width=True,
            on_select="rerun", selection_mode="points", key=f"pd_chart_{liquid}")

        # ── Transition jump (ΔV / ΔS / ΔH) on click ───────────────────────
        _render_transition_jump(event)

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
            st.download_button("Download CSV", make_csv(export_df),
                               file_name="seafreeze_phase_diagram.csv",
                               mime="text/csv")
