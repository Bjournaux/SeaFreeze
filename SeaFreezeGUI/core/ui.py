"""Shared UI helpers for the SeaFreeze GUI views."""

import os
import sys

import pandas as pd
import plotly.graph_objects as go

from core.constants import MATERIAL_LABELS, MATERIAL_SHORT_LABELS

# ─────────────────────────────────────────────────────────────────────────────
# Asset path resolution (works for both `streamlit run` and PyInstaller bundle)
# ─────────────────────────────────────────────────────────────────────────────
_ASSETS_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "assets")

def asset(name):
    if hasattr(sys, "_MEIPASS"):
        return os.path.join(sys._MEIPASS, "assets", name)
    return os.path.join(_ASSETS_DIR, name)


# ─────────────────────────────────────────────────────────────────────────────
# Labels & formatting
# ─────────────────────────────────────────────────────────────────────────────

def label(material):
    """Full citation label — used in selectors and page headers."""
    return MATERIAL_LABELS.get(material, material)

def short_label(material):
    """Short label — used in plot legends to keep trace names compact."""
    return MATERIAL_SHORT_LABELS.get(material, material)

def fmt_range(lo, hi, unit):
    return f"{lo:.1f} – {hi:.1f} {unit}"

def make_csv(df: pd.DataFrame) -> str:
    return df.to_csv(index=False)

def title_with_fixed(name, const_vars):
    """Build plot title appending fixed variable values."""
    if not const_vars:
        return name
    fixed = ", ".join(f"{v[0]} = {v[1]:.4g} {v[2]}" for v in const_vars)
    return f"{name}<br><sup>{fixed}</sup>"


# ─────────────────────────────────────────────────────────────────────────────
# Stability-field boundary overlays
# ─────────────────────────────────────────────────────────────────────────────
BOUNDARY_COLORS = [
    "#FF6347", "#FFD700", "#00FFFF", "#FF69B4",
    "#7FFF00", "#FF8C00", "#DA70D6", "#1E90FF",
]


def add_phase_boundaries(fig, phase_bounds, var1_name, var2_name, mode="2d"):
    """Overlay phase boundary lines on a Plotly figure."""
    for i, (mat, other, bP, bT) in enumerate(phase_bounds):
        if var1_name == "P":
            x_vals, y_vals = bP, bT
        else:
            x_vals, y_vals = bT, bP
        color = BOUNDARY_COLORS[i % len(BOUNDARY_COLORS)]
        name = f"{mat}–{other}"
        if mode == "3d":
            fig.add_trace(go.Scatter3d(
                x=x_vals, y=y_vals,
                z=[None] * len(x_vals),
                mode="lines",
                line=dict(color=color, width=4),
                name=name,
                showlegend=True,
            ))
        else:
            fig.add_trace(go.Scatter(
                x=x_vals, y=y_vals,
                mode="lines",
                line=dict(color=color, width=2.5, dash="dot"),
                name=name,
                showlegend=True,
            ))
