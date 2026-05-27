# SeaFreeze GUI

Interactive web and desktop application for computing thermodynamic properties of water, ice polymorphs, and NaCl aqueous solutions using [SeaFreeze](https://github.com/Bjournaux/SeaFreeze).

Built with [Streamlit](https://streamlit.io/) and [Plotly](https://plotly.com/python/).

---

## Features

### Property Calculator

Compute thermodynamic properties for any SeaFreeze material (Ice Ih–VI, Water, NaCl(aq)) with adaptive visualization that adjusts to your input:

- **Single point** — returns a table of all properties with values and units
- **1-D sweep** — vary one coordinate (P, T, or m) to produce interactive line plots
- **2-D grid** — vary two coordinates to produce heatmaps, heatmaps with isocontours, or 3-D surface plots

All modes support CSV export. The 2-D mode includes selectable color scales, adjustable contour density, and optional stability field boundary overlays.

| Single point | 1-D sweep |
|:---:|:---:|
| ![Single point](screenshots/Ice%20VI%20Single%20PT.png) | ![1-D sweep](screenshots/Ice%20Ih%20T%20constant%2C%20P%20range.png) |

| Heatmap | Heatmap + Isocontours | 3-D Surface |
|:---:|:---:|:---:|
| ![Heatmap](screenshots/Ice%20V%2C%20P%20and%20T%20range%20Heatmap.png) | ![Isocontours](screenshots/Ice%20V%2C%20P%20and%20T%20range%20Heatmap%2Bisocontour.png) | ![3D](screenshots/Ice%20V%2C%20P%20and%20T%20range%203D%20surface.png) |

### Phase Diagram

Interactive water phase diagram with selectable phase boundaries:

- Checkboxes for each ice polymorph (Ih, II, III, V, VI) and liquid water
- Toggle between stable, metastable, or all segments
- Overlay NaCl(aq) melting curves at user-specified molalities
- Adjustable P and T axis ranges
- CSV export of all displayed boundary data

| Pure water phase diagram | With NaCl(aq) melting curves |
|:---:|:---:|
| ![Phase diagram](screenshots/Phase%20diagram.png) | ![Phase diagram + NaCl](screenshots/Phase%20diagram%20%2B%202m%20NaCl.png) |

---

## Quick Start (Web App)

### Requirements

- Python 3.10+
- pip

### Install and run

```bash
# Clone the repository
git clone https://github.com/Bjournaux/SeaFreeze.git
cd SeaFreeze/SeaFreezeGUI

# Install dependencies
pip install -r requirements.txt

# Launch the app
streamlit run app.py
```

The app opens in your browser at `http://localhost:8501`.

### Dependencies

| Package | Version |
|---------|---------|
| SeaFreeze | >= 1.1.1 |
| Streamlit | >= 1.39 |
| Plotly | >= 5.20 |
| Pandas | >= 2.1 |
| NumPy | >= 1.26 |

---

## Building the Desktop App (macOS Apple Silicon)

The desktop app bundles the full Python interpreter and all dependencies into a standalone `.app` — no Python installation needed on the target machine.

### Prerequisites

```bash
pip install pyinstaller
```

### Build

```bash
cd SeaFreezeGUI
python -m PyInstaller SeaFreezeGUI.spec --noconfirm
```

The output is at `dist/SeaFreezeGUI.app` (~335 MB). Double-click to launch — it starts a local Streamlit server and opens your browser automatically.

### Notes

- **First launch on macOS:** The app is unsigned, so Gatekeeper will block it. Right-click the `.app` and select "Open" to bypass, or sign it with an Apple Developer certificate.
- **Rebuilding after changes:** Delete `build/` and `dist/`, then re-run the PyInstaller command.
- The build spec (`SeaFreezeGUI.spec`) is pre-configured for ARM64 (Apple Silicon). For Intel Macs, change `target_arch` to `"x86_64"`.

---

## Project Structure

```
SeaFreezeGUI/
├── app.py                  # Main Streamlit app (Property Calculator + Phase Diagram)
├── launcher.py             # Desktop launcher (starts Streamlit in-process)
├── SeaFreezeGUI.spec       # PyInstaller build specification
├── requirements.txt        # Python dependencies
├── core/
│   ├── __init__.py
│   ├── compute.py          # Cached wrappers around SeaFreeze API
│   └── constants.py        # Material lists, property metadata, categories
├── .streamlit/
│   └── config.toml         # Streamlit theme and server config
├── screenshots/            # App screenshots for documentation
└── assets/                 # Icons and static assets
```

---

## Citation

If you use SeaFreeze in your research, please cite:

Journaux et al. (2020). Holistic Approach for Studying Planetary Hydrospheres: Gibbs Representation of Ices Thermodynamics, Elasticity, and the Water Phase Diagram to 2300 MPa. *Journal of Geophysical Research: Planets*, 125, e2019JE006176. https://doi.org/10.1029/2019JE006176
