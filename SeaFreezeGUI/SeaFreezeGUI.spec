# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec for SeaFreeze GUI — Apple Silicon (arm64).

Bundles the full Python interpreter + all dependencies.
No external Python installation needed on the target machine.
"""

import os
import importlib
from PyInstaller.utils.hooks import collect_all, collect_data_files, collect_submodules

block_cipher = None


# ── Helper ────────────────────────────────────────────────────────────────────
def pkg_path(name):
    mod = importlib.import_module(name)
    return os.path.dirname(mod.__file__)


# ── Collect full packages (data + submodules) ─────────────────────────────────
# Streamlit has many runtime data files (static/, proto/, etc.)
st_datas, st_binaries, st_hiddenimports = collect_all("streamlit")

# Plotly templates and data
plotly_datas = collect_data_files("plotly")

# SeaFreeze spline data
sf_datas = collect_data_files("seafreeze")

# ── Assemble datas ────────────────────────────────────────────────────────────
datas = [
    # App files
    ("app.py", "."),
    ("core", "core"),
    (".streamlit", ".streamlit"),
]

# Add assets if they exist
if os.path.isdir("assets"):
    datas.append(("assets", "assets"))

datas += st_datas + plotly_datas + sf_datas

# ── Hidden imports ────────────────────────────────────────────────────────────
hiddenimports = [
    # Streamlit internals (it uses dynamic imports heavily)
    *st_hiddenimports,
    *collect_submodules("streamlit"),
    *collect_submodules("plotly"),

    # SeaFreeze and its dependencies
    "seafreeze",
    "seafreeze.seafreeze",
    "seafreeze.phaselines",
    "lbftd",
    "lbftd.fdgrid",
    "mlbspline",
    "mlbspline.load",

    # Scientific stack
    "h5py",
    "scipy",
    "scipy.interpolate",
    "scipy.special",
    "scipy.special._cdflib",
    "numpy",
    "numpy.core",
    "pandas",
    "pandas._libs",

    # Streamlit runtime deps
    "pyarrow",
    "pyarrow.lib",
    "PIL",
    "toml",
    "click",
    "altair",
    "jinja2",
    "markupsafe",
    "packaging",
    "importlib_metadata",
    "pydeck",
    "tornado",
    "tornado.web",
    "tornado.websocket",
    "watchdog",
    "cachetools",
    "validators",
    "gitdb",
    "tenacity",
    "rich",
    "typing_extensions",
]

# ── Analysis ──────────────────────────────────────────────────────────────────
a = Analysis(
    ["launcher.py"],
    pathex=[],
    binaries=st_binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        "tkinter",
        "matplotlib",
        "IPython",
        "notebook",
        "pytest",
    ],
    noarchive=False,
    optimize=0,
    cipher=block_cipher,
)

pyz = PYZ(a.pure, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="SeaFreezeGUI",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,       # No terminal window
    target_arch="arm64",
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    name="SeaFreezeGUI",
)

app = BUNDLE(
    coll,
    name="SeaFreezeGUI.app",
    icon=None,  # Set to "assets/icon.icns" for a custom icon
    bundle_identifier="edu.uw.seafreeze.gui",
    info_plist={
        "CFBundleShortVersionString": "1.1.2",
        "CFBundleName": "SeaFreeze GUI",
        "NSHighResolutionCapable": True,
    },
)
