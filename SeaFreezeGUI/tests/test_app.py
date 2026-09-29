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
@pytest.mark.parametrize("material", ["Ih", "water1", "NaClaq"])
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
    assert {"Ih–water1", "Ih–II", "Ih–III"} <= set(names)


def test_three_ranges_warns():
    at = _open()
    at.selectbox(key="pc_material").set_value("NaClaq")
    _run(at)
    for var in ("P", "T", "m"):
        _set_range(at, var, 5)
    _run(at)
    assert any("3 independent ranges" in w.value for w in at.warning)


# ── Phase Diagram ────────────────────────────────────────────────────────────
def test_phase_diagram_default():
    at = _open("Phase Diagram")
    [fig] = _figures(at)
    names = _trace_names(fig)
    assert "Ice Ih – Water" in names
    assert "Ice VI – Water" in names
    assert fig["layout"]["xaxis"]["title"]["text"] == "Pressure (MPa)"
    texts = [a["text"] for a in fig["layout"].get("annotations", [])]
    assert "<b>Liquid</b>" in texts
    labels = [b.label for b in at.get("download_button")]
    assert "Download CSV" in labels


def test_phase_diagram_metastable_and_triple_points():
    at = _open("Phase Diagram")
    at.radio(key="pd_segment").set_value("all")
    at.checkbox(key="pd_show_tp").check()
    _run(at)
    [fig] = _figures(at)
    names = _trace_names(fig)
    assert "Triple points" in names
    assert any(n.endswith("(meta)") for n in names)


def test_phase_diagram_nacl_curves():
    at = _open("Phase Diagram")
    at.checkbox(key="pd_show_nacl").check()
    _run(at)
    assert at.text_input(key="pd_nacl_m").value == "1, 2, 3"
    [fig] = _figures(at)
    names = _trace_names(fig)
    for m in (1.0, 2.0, 3.0):
        assert f"Ice Ih – NaCl(aq) m={m}" in names


def test_phase_diagram_too_few_phases():
    at = _open("Phase Diagram")
    for phase in ["II", "III", "V", "VI", "water1"]:
        at.checkbox(key=f"pd_{phase}").uncheck()
    _run(at)
    assert not _figures(at)
    assert any("at least two phases" in i.value for i in at.info)


# ── About ────────────────────────────────────────────────────────────────────
def test_about_table():
    at = _open("About")
    assert len(at.dataframe) == 1
    df = at.dataframe[0].value
    assert list(df["Material"]) == [MATERIAL_LABELS[m] for m in ALL_MATERIALS]
    assert {"P range (MPa)", "T range (K)"} <= set(df.columns)
