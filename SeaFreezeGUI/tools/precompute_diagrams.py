"""Precompute the default full phase diagrams shipped with the GUI.

Writes assets/diagrams/phase_diagram_PT.npz and phase_diagram_rhoT.npz (the
default window at the default resolution, see core/diagrams.py).  The GUI
checks the stored window, resolution and SeaFreeze version, and recomputes
live if they do not match — so rerun this after changing any of them:

    cd SeaFreezeGUI && python3 tools/precompute_diagrams.py
"""
import os
import sys
import time

GUI = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, GUI)
sys.path.insert(0, os.path.join(os.path.dirname(GUI), "Python"))

from core import diagrams as D  # noqa: E402

for coords in ("PT", "rhoT"):
    window = D.DEFAULT_WINDOW[coords]
    t = time.time()
    d = D.build_diagram(coords, window, D.DEFAULT_RESOLUTION)
    path = D.precomputed_path(coords)
    D.save_diagram(d, path, coords, window, D.DEFAULT_RESOLUTION)
    back = D.load_diagram(path, coords, window, D.DEFAULT_RESOLUTION)
    assert back is not None, "saved diagram does not load back"
    print(f"{coords:5s} {time.time() - t:5.1f} s  {os.path.getsize(path) / 1e6:.2f} MB  -> {path}")
