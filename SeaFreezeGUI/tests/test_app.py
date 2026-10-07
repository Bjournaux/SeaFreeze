"""Smoke tests for the SeaFreeze GUI, driven by streamlit.testing.v1.AppTest.

Run from the GUI folder:  cd SeaFreezeGUI && python3 -m pytest -q tests

The app finds the in-repo SeaFreeze library (../Python) on its own; outside a
repo checkout, SeaFreeze must be installed.
"""

import json
import os
import sys

import pytest
from streamlit.testing.v1 import AppTest

GUI_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
APP_PATH = os.path.join(GUI_DIR, "app.py")
if GUI_DIR not in sys.path:
    sys.path.insert(0, GUI_DIR)

from core.constants import ALL_MATERIALS, MATERIAL_LABELS  # noqa: E402

# First runs load splines and contour phase lines — allow plenty of time.
TIMEOUT = 120

PAGES = ["Property Calculator", "Phase Diagram", "About"]


# ── Helpers ──────────────────────────────────────────────────────────────────
def _run(at):
    at.run(timeout=TIMEOUT)
    assert not at.exception, [e.value for e in at.exception]
    return at


def _open(page="Property Calculator"):
    at = _run(AppTest.from_file(APP_PATH, default_timeout=TIMEOUT))
    if page != "Property Calculator":
        at.pills(key="page_select").set_value(page)
        _run(at)
    return at


def _melting(at):
    """Switch the Phase Diagram page to the ice melting-line viewer."""
    at.radio(key="pd_mode").set_value("Ice melting lines")
    return _run(at)


def _compute(at):
    [btn] = [b for b in at.button if b.label == "Compute"]
    btn.click()
    return _run(at)


def _figures(at):
    """Plotly figure specs (dicts with 'data' and 'layout') on the page."""
    return [json.loads(c.proto.spec) for c in at.get("plotly_chart")]


def _trace_names(fig):
    return [tr.get("name") for tr in fig["data"]]


def _set_range(at, var, npts):
    at.radio(key=f"{var}_mode").set_value("Range")
    _run(at)
    at.number_input(key=f"{var}_npts").set_value(npts)


# ── Every page renders ───────────────────────────────────────────────────────
@pytest.mark.parametrize("page", PAGES)
def test_page_renders(page):
    at = _open(page)
    assert at.title[0].value.startswith("SeaFreeze")
    assert at.pills(key="page_select").value == page


# ── Property Calculator ──────────────────────────────────────────────────────
@pytest.mark.parametrize("material", ["Ih", "water_Bollengier2019", "NaClaq_Brown2026"])
def test_single_point(material):
    at = _open()
    at.selectbox(key="pc_material").set_value(material)
    _run(at)
    assert "Single point" in at.info[0].value
    _compute(at)

    assert len(at.dataframe) == 1
    df = at.dataframe[0].value
    assert len(df) > 0
    assert list(df.columns) == ["Property", "Description", "Value", "Unit"]
    assert "rho" in set(df["Property"])
    assert at.subheader[0].value == f"Results — {MATERIAL_LABELS[material]}"


def test_single_point_selected_properties():
    at = _open()
    at.checkbox(key="select_all").uncheck()
    _run(at)
    assert at.multiselect(key="prop_select").value == ["rho", "Cp", "G"]
    _compute(at)
    df = at.dataframe[0].value
    assert set(df["Property"]) == {"rho", "Cp", "G"}


def test_1d_sweep():
    at = _open()
    at.selectbox(key="pc_material").set_value("Ih")
    _set_range(at, "P", 12)
    _run(at)
    assert "1-D sweep" in at.info[0].value
    _compute(at)

    figs = _figures(at)
    assert len(figs) > 1                       # one line plot per property tab
    for fig in figs:
        assert fig["layout"]["xaxis"]["title"]["text"] == "P (MPa)"
        assert len(fig["data"][0]["x"]) == 12
    labels = [b.label for b in at.get("download_button")]
    assert "Download CSV" in labels


@pytest.mark.parametrize("view_mode",
                         ["Heatmap", "Heatmap + Isocontours", "3-D Surface"])
def test_2d_grid_with_boundaries(view_mode):
    at = _open()
    at.selectbox(key="pc_material").set_value("Ih")
    _set_range(at, "P", 15)
    _set_range(at, "T", 15)
    _run(at)
    assert "2-D grid" in at.info[0].value
    _compute(at)

    assert at.selectbox(key="prop_2d").value is not None
    at.radio(key="view_mode").set_value(view_mode)
    at.checkbox(key="show_boundaries").check()
    _run(at)

    [fig] = _figures(at)
    types = [tr["type"] for tr in fig["data"]]
    names = _trace_names(fig)
    main = {"Heatmap": "heatmap", "Heatmap + Isocontours": "heatmap",
            "3-D Surface": "surface"}[view_mode]
    assert types[0] == main
    if view_mode == "Heatmap + Isocontours":
        assert "contour" in types
    # Ice Ih borders liquid water, ice II and ice III
    assert {"Ih–water_Bollengier2019", "Ih–II", "Ih–III"} <= set(names)


def test_three_ranges_warns():
    at = _open()
    at.selectbox(key="pc_material").set_value("NaClaq_Brown2026")
    _run(at)
    for var in ("P", "T", "m"):
        _set_range(at, var, 5)
    _run(at)
    assert any("3 independent ranges" in w.value for w in at.warning)


# ── Phase Diagram ────────────────────────────────────────────────────────────
def test_phase_diagram_default():
    at = _melting(_open("Phase Diagram"))
    [fig] = _figures(at)
    names = _trace_names(fig)
    assert "Ice Ih – Water (Bollengier 2019)" in names
    assert "Ice VI – Water (Bollengier 2019)" in names
    assert fig["layout"]["xaxis"]["title"]["text"] == "Pressure (MPa)"
    texts = [a["text"] for a in fig["layout"].get("annotations", [])]
    assert "<b>Liquid</b>" in texts
    labels = [b.label for b in at.get("download_button")]
    assert "Download CSV" in labels


def test_phase_diagram_metastable_and_triple_points():
    at = _melting(_open("Phase Diagram"))
    at.radio(key="pd_segment").set_value("all")
    at.checkbox(key="pd_show_tp").check()
    _run(at)
    [fig] = _figures(at)
    names = _trace_names(fig)
    assert "Triple points" in names
    assert any(n.endswith("(meta)") for n in names)


def test_phase_diagram_nacl_curves():
    at = _melting(_open("Phase Diagram"))
    at.checkbox(key="pd_show_nacl").check()
    _run(at)
    assert at.text_input(key="pd_nacl_m").value == "1, 2, 3"
    [fig] = _figures(at)
    names = _trace_names(fig)
    for m in (1.0, 2.0, 3.0):
        assert f"Ice Ih – NaCl(aq) m={m}" in names


def test_phase_diagram_too_few_phases():
    at = _melting(_open("Phase Diagram"))
    for key in ["pd_II", "pd_III", "pd_V", "pd_VI", "pd_liquid_on"]:
        at.checkbox(key=key).uncheck()
    _run(at)
    assert not _figures(at)
    assert any("at least two phases" in i.value for i in at.info)


def test_melting_lines_water_Brown2026_and_log_axes():
    at = _melting(_open("Phase Diagram"))
    at.radio(key="pd_liquid").set_value("water_Brown2026")
    at.checkbox(key="pd_logP").check()
    at.checkbox(key="pd_logT").check()
    _run(at)
    [fig] = _figures(at)
    names = _trace_names(fig)
    assert "Ice Ih – Water (Brown 2026)" in names and "Ice VI – Water (Brown 2026)" in names
    assert fig["layout"]["xaxis"]["type"] == "log" and fig["layout"]["yaxis"]["type"] == "log"


# ── Full diagram ─────────────────────────────────────────────────────────────
def test_full_diagram_default_is_precomputed():
    at = _open("Phase Diagram")
    assert at.radio(key="pd_mode").value.startswith("Full diagram")
    assert any("precomputed default" in c.value for c in at.caption)
    [fig] = _figures(at)
    names = _trace_names(fig)
    assert fig["data"][0]["type"] == "heatmap"
    assert {"phase boundaries", "saturation curve", "critical point", "triple points"} <= set(names)
    assert fig["layout"]["xaxis"]["type"] == "log"
    texts = {a["text"] for a in fig["layout"]["annotations"]}
    assert {"<b>liquid</b>", "<b>vapour</b>", "<b>VI</b>"} <= texts


@pytest.mark.parametrize("coords", ["P–T", "ρ–T"])
@pytest.mark.parametrize("view", ["2-D map", "2-D map + isocontours", "3-D surface"])
def test_full_diagram_property_views(coords, view):
    at = _open("Phase Diagram")
    at.radio(key="pdf_coords").set_value(coords)
    _run(at)
    at.selectbox(key="pdf_colour").set_value("rho")
    _run(at)
    at.radio(key="pdf_view").set_value(view)
    at.checkbox(key="pdf_ylog").check()
    _run(at)
    [fig] = _figures(at)
    types = [tr["type"] for tr in fig["data"]]
    if view == "3-D surface":
        assert types[0] == "surface"
        assert fig["layout"]["scene"]["yaxis"]["type"] == "log"
    else:
        assert "heatmap" in types and fig["layout"]["yaxis"]["type"] == "log"
        assert ("contour" in types) == (view == "2-D map + isocontours")
    if coords == "ρ–T" and view != "3-D surface":
        assert "saturation dome" in _trace_names(fig)


def test_full_diagram_linear_axes():
    at = _open("Phase Diagram")
    at.checkbox(key="pdf_xlog_PT").uncheck()
    _run(at)
    [fig] = _figures(at)
    assert fig["layout"]["xaxis"]["type"] == "linear"


# ── water_Brown2026 and (rho, T) input in the calculator ─────────────────────────────
@pytest.mark.parametrize("branch, lo, hi", [("stable", 990, 1000), ("liquid", 990, 1000),
                                            ("vapor", 0, 0.1)])
def test_water_Brown2026_single_point_branches(branch, lo, hi):
    at = _open()
    at.selectbox(key="pc_material").set_value("water_Brown2026")
    _run(at)
    at.number_input(key="P_single_water_Brown2026_PT").set_value(0.001)
    at.number_input(key="T_single_water_Brown2026_PT").set_value(300.0)
    at.radio(key="pc_branch").set_value(branch)
    _run(at)
    _compute(at)
    df = at.dataframe[0].value
    rho = float(df.set_index("Property").loc["rho", "Value"])
    if branch == "stable":                  # 1 kPa < p_sat(300 K): the vapour is stable
        lo, hi = 0, 0.1
    assert lo < rho < hi


def test_water_Brown2026_rhoT_log_sweep():
    at = _open()
    at.selectbox(key="pc_material").set_value("water_Brown2026")
    _run(at)
    at.radio(key="pc_input").set_value("ρ, T")
    _run(at)
    _set_range(at, "rho", 20)
    _run(at)
    assert at.checkbox(key="rho_log_water_Brown2026_rhoT").value     # log by default for water_Brown2026
    _compute(at)
    figs = _figures(at)
    assert figs and all(f["layout"]["xaxis"]["type"] == "log" for f in figs)
    assert figs[0]["layout"]["xaxis"]["title"]["text"] == "ρ (kg/m³)"
    # pressure is an output in (rho, T) mode
    assert "P" in [t.label for t in at.tabs]


def test_water_Brown2026_2d_with_phase_diagram_overlay():
    at = _open()
    at.selectbox(key="pc_material").set_value("water_Brown2026")
    _run(at)
    _set_range(at, "P", 30)
    _set_range(at, "T", 20)
    _run(at)
    _compute(at)
    at.checkbox(key="show_boundaries").check()
    _run(at)
    [fig] = _figures(at)
    names = _trace_names(fig)
    assert fig["layout"]["xaxis"]["type"] == "log"
    assert {"saturation curve", "critical point"} <= set(names)


def test_gibbs_rhoT_input_small_grid():
    at = _open()
    at.selectbox(key="pc_material").set_value("Ih")
    _run(at)
    at.radio(key="pc_input").set_value("ρ, T")
    _run(at)
    at.number_input(key="rho_single_Ih_rhoT").set_value(920.0)
    at.number_input(key="T_single_Ih_rhoT").set_value(260.0)
    _run(at)
    _compute(at)
    df = at.dataframe[0].value.set_index("Property")
    P = float(df.loc["P", "Value"])
    assert 0 < P < 100                       # ice Ih at 920 kg/m3, 260 K


# ── About ────────────────────────────────────────────────────────────────────
def test_about_table():
    at = _open("About")
    assert len(at.dataframe) == 1
    df = at.dataframe[0].value
    assert list(df["Material"]) == [MATERIAL_LABELS[m] for m in ALL_MATERIALS]
    assert {"P range (MPa)", "T range (K)"} <= set(df.columns)
