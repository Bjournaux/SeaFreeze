"""Test setup for the SeaFreeze GUI.

Run from the repository root:

    python -m pytest SeaFreezeGUI/tests

Puts the GUI directory and the in-repo Python package on sys.path so the
tests exercise the checked-out seafreeze code, not an installed copy.
"""

import os
import sys

import pytest

_GUI_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_PY_DIR = os.path.join(os.path.dirname(_GUI_DIR), 'Python')
for p in (_PY_DIR, _GUI_DIR):
    if p not in sys.path:
        sys.path.insert(0, p)

streamlit = pytest.importorskip('streamlit', reason='GUI tests need streamlit installed')
