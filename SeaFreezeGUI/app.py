"""SeaFreeze GUI — Streamlit application.

Property Calculator + Phase Diagram viewer. Each page lives in views/.
"""

import sys, os
# Ensure core/ and views/ are importable — works from source and PyInstaller bundle
_app_dir = os.path.dirname(os.path.abspath(__file__))
if _app_dir not in sys.path:
    sys.path.insert(0, _app_dir)
# PyInstaller bundles files under sys._MEIPASS
if hasattr(sys, "_MEIPASS") and sys._MEIPASS not in sys.path:
    sys.path.insert(0, sys._MEIPASS)
# From a repo checkout, use the in-repo SeaFreeze library (<repo>/Python)
# ahead of any installed copy; otherwise fall back to the installed package.
_repo_python = os.path.join(os.path.dirname(_app_dir), "Python")
if (not hasattr(sys, "_MEIPASS")
        and os.path.isfile(os.path.join(_repo_python, "seafreeze", "__init__.py"))
        and _repo_python not in sys.path):
    sys.path.insert(0, _repo_python)

import streamlit as st

from core.ui import asset
from views import about, phase_diagram, property_calculator

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
# Sidebar branding — shown on every page
# ─────────────────────────────────────────────────────────────────────────────
with st.sidebar:
    logo_path = asset("logo.png")
    if os.path.exists(logo_path):
        import base64 as _b64
        with open(logo_path, "rb") as _f:
            _logo_b64 = _b64.b64encode(_f.read()).decode()
        st.markdown(
            f'<a href="https://bjournaux.wordpress.com/" target="_blank">'
            f'<img src="data:image/png;base64,{_logo_b64}" width="230"></a>',
            unsafe_allow_html=True,
        )
    st.caption(
        "Developed by **Baptiste Journaux**  \n"
        "University of Washington"
    )
    st.divider()


# ─────────────────────────────────────────────────────────────────────────────
# Page selector (top of main area) and routing
# ─────────────────────────────────────────────────────────────────────────────
_PAGES = {
    "Property Calculator": property_calculator,
    "Phase Diagram":       phase_diagram,
    "About":               about,
}

page = st.pills("Navigation",
                list(_PAGES),
                default="Property Calculator",
                key="page_select", label_visibility="collapsed")

# A deselected pill (None) falls through to About
_PAGES.get(page, about).render()
