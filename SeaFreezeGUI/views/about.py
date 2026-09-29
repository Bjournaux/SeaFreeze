"""About page — validity ranges, validation figure, contributors and references."""

import os

import pandas as pd
import streamlit as st

from core.constants import ALL_MATERIALS, MATERIAL_LABELS
from core.compute import get_phase_range
from core.ui import asset


def render():
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
    fig_path = asset("Phase_diagram_exp_data.png")
    if os.path.exists(fig_path):
        st.image(fig_path,
                 caption="SeaFreeze computed phase boundaries vs. experimental data.",
                 use_container_width=True)
    else:
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
