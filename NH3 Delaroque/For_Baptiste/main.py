####################
# Import  librairies
####################

import numpy as np
import matplotlib.pyplot as plt

from utils import *

from scipy.optimize import least_squares

####################
#       Ice Ih
#####################
# Here, fitting data and find L_i

interp_liq, interp_ih = initialize_fusion_model(
    "data/seafreeze/GibbsEnergywater_IAPWS95_kg.txt",
    "data/seafreeze/GibbsEnergyIh_kg.txt"
)

# Mass fractions are converted into molar fractions
datasets = [
    ("data/literature/amm2pct_ih.csv", 0.0211),
    ("data/literature/amm8pct_ih.csv", 0.0842),
    ("data/literature/amm11pct_ih.csv", 0.1156),
    ("data/literature/Leliwa_x3.34.csv", 0.0334),
    ("data/literature/Leliwa_x4.72.csv", 0.0472),
    ("data/literature/Leliwa_x11.1.csv", 0.111),
    ("data/literature/Leliwa_x17.6.csv", 0.1843),
]

P_all, T_all, x_all = load_liquidus_datasets(datasets)

# Fitting of the Redlich-Kister parameters
fit = fit_margules(P_all, T_all, x_all, interp_liq, interp_ih)

L_fit = fit.x

if not fit.success:
    raise RuntimeError(f"Fit for ice Ih failed: {fit.message}")

print("Fitted Relich-Kister parameters:")
print(f"L0 = {L_fit[0]:.6f}")
print(f"L1 = {L_fit[1]:.6f}")
print(f"L2 = {L_fit[2]:.6f}")

res = residuals(fit.x, P_all, T_all, x_all, interp_liq, interp_ih)


##############################
#            Ice V
##############################

interp_liq_V, interp_V = initialize_fusion_model(
    "data/seafreeze/GibbsEnergywater_IAPWS95_kg_V.txt",
    "data/seafreeze/GibbsEnergyV_kg.txt"
)

datasets_V = [
    ("data/literature/amm015_V.csv", 0.0741),
    ("data/literature/Chroukroun_x01_V.csv", 0.1052),
    ("data/literature/Chroukroun_x02_V.csv", 0.209),
]

PV_all, TV_all, xV_all = load_liquidus_datasets(datasets_V)

fit_V = least_squares(residuals_V, x0=[0.0],
    args=(PV_all, TV_all, xV_all, interp_liq_V, interp_V, L_fit))

if not fit_V.success:
    raise RuntimeError(f"Fit for ice V failed: {fit_V.message}")

C0_fit = fit_V.x[0]

print(f"C0 (fit V) = {C0_fit:.6f}")



##############################
#            Ice III
##############################

# Initialize Gibbs energy interpolators for the liquid and ice III phases.
interp_liq_III, interp_III = initialize_fusion_model(
    "data/seafreeze/GibbsEnergywater_IAPWS95_kg_III.txt",
    "data/seafreeze/GibbsEnergyIII_kg.txt"
)

datasets_III = [
    ("data/literature/amm014_III.csv", 0.0264),
    ("data/literature/amm015_III.csv", 0.0741),
    ("data/literature/Chroukroun_x01_III.csv", 0.1052),
]

PIII_all, TIII_all, xIII_all = load_liquidus_datasets(datasets_III)

# Pressure range [MPa]
PIII_model = np.linspace(215, 378, 80)

# NH3 compositions of the data
x_values = [0.1052, 0.0741, 0.209]

# Calculate liquidus curves
TIII_curves = {}

for x in x_values:

    TIII_curves[x] = np.array([find_liquidus_temperature(P, x, interp_liq_III, interp_III,  
            L_fit=L_fit, T_bounds=(180, 300)) for P in PIII_model])

TIII_x01  = TIII_curves[0.1052]
TIII_x073 = TIII_curves[0.0741]
TIII_x02  = TIII_curves[0.209]
mask = np.isfinite(TIII_x01)


##############################
#            Ice II
##############################

interp_liq_II, interp_II = initialize_fusion_model(
    "data/seafreeze/GibbsEnergywater_IAPWS95_kg_II.txt",
    "data/seafreeze/GibbsEnergyII_kg.txt"
)

PII_model = np.linspace(192, 503, 100)

TII_curves = {}

for x in x_values:

    TII_curves[x] = np.array([find_liquidus_temperature(P, x, interp_liq_II, interp_II,
            L_fit=L_fit, T_bounds=(180, 300)) for P in PII_model])

TII_x01  = TII_curves[0.1052]
TII_x073 = TII_curves[0.0741]
TII_x02  = TII_curves[0.209]

##############################
#            Ice VI
##############################

interp_liq_VI, interp_VI = initialize_fusion_model(
    "data/seafreeze/GibbsEnergywater_IAPWS95_kg_VI.txt",
    "data/seafreeze/GibbsEnergyVI_kg.txt"
)

PVI_model = np.linspace(650, 900, 100)

TVI_curves = {}

for x in x_values:

    TVI_curves[x] = np.array([find_liquidus_temperature(P, x, interp_liq_VI, interp_VI,
            L_fit=L_fit, T_bounds=(180, 300)) for P in PVI_model])

TVI_x01  = TVI_curves[0.1052]
TVI_x073 = TVI_curves[0.0741]
TVI_x02  = TVI_curves[0.209]

##############################
#       Ice models
##############################

# initialization
ice_models = {
    "Ih": {
        "liquid": interp_liq,
        "solid": interp_ih,
        "C0": 0.0,
        "C1": 0.0,
        "C2": 0.0,
    },
    "II": {
        "liquid": interp_liq_II,
        "solid": interp_II,
        "C0": 0.0,
        "C1": 0.0,
        "C2": 0.0,
    },
    "III": {
        "liquid": interp_liq_III,
        "solid": interp_III,
        "C0": 0.0,
        "C1": 0.0,
        "C2": 0.0,
    },
    "V": {
        "liquid": interp_liq_V,
        "solid": interp_V,
        "C0": 0.0,  
        "C1": 0.0,
        "C2": 0.0,
    },
    "VI": {
        "liquid": interp_liq_VI,
        "solid": interp_VI,
        "C0": 0.0,
        "C1": 0.0,
        "C2": 0.0,
    },
}


########################################
# Determine Tliq with phase stability
########################################

Pgrid = np.linspace(0, 900, 400)

compositions = [0.1052, 0.209]

results = {}

for x in compositions:

    results[x] = [liquidus_global(P, x, ice_models, L_fit) for P in Pgrid]

T_curves = {x: np.array([result[0] for result in results[x]]) for x in compositions}
phase_curves = {x: np.array([result[1] for result in results[x]]) for x in compositions}


########################################
# Stable ice phase from pure-water map
########################################

transitions = ["Ih_II", "Ih_III", "II_III",
    "II_V", "II_VI", "III_V", "V_VI",]

phase_boundaries = {}
triple_points = {}

for transition in transitions:

    T_boundary, P_boundary = load_phase_boundary(f"data/seafreeze/{transition}.txt")
    triples = load_triple_points(f"data/seafreeze/TriplePoints_{transition}.txt")

    T_boundary, P_boundary = crop_transition(T_boundary, P_boundary, triples)
    phase_boundaries[transition] = {"T": T_boundary, "P": P_boundary,}

    triple_points[transition] = triples

########################################
#         Final phase diagram
########################################

Pgrid = np.linspace(0, 900, 400)
compositions = [0.0211, 0.0842, 0.1156, 0.1843]

# Pure water
pure_water_results = [liquidus_global(P, 0.0, ice_models, L_fit) for P in Pgrid]

T_pure = np.array([result[0] for result in pure_water_results])

# storing the curves of several composition of the H2O-NH3 mixture
T_curves = {}
phase_curves = {}

for x in compositions:

    results = [liquidus_global(P, x, ice_models, L_fit) for P in Pgrid]

    T_curves[x] = np.array([result[0] for result in results])
    phase_curves[x] = np.array([result[1] for result in results])


# Plotting
fig, ax = plt.subplots()

ax.plot(Pgrid, T_pure, c="black", lw=2, label="Pure water")

for x in compositions:

    ax.plot(Pgrid,  T_curves[x], label=f"x(NH$_3$) = {100*x:.0f}%")

# Pure-water solid-solid phase boundaries
for transition, boundary in phase_boundaries.items():

    ax.plot(boundary["P"], boundary["T"], c="gray", lw=1, alpha=0.7, ls=":",
        label="Pure-water phase boundaries" if transition == transitions[0] else None)

ax.set_xlabel("Pressure [MPa]", fontsize=14)
ax.set_ylabel("Temperature [K]", fontsize=14)

ax.set_xlim(0, 1000)
ax.set_ylim(200, 300)

# Datapoints for comparison with model prediction!!!
# PLEASE COMMENT OR UNCOMMENT THOSE LINES

# comp_data = np.array([x for _, x in datasets])
# colors = ['dodgerblue', 'green', 'lightgray',
#           'lightgray', 'lightgray', 'lightgray', 'r']

# markers = ['o', 's', '^', 's', '+', '+', '+']

# data_lengths = [len(load_liquidus_datasets([d])[0]) for d in datasets]

# start = 0

# for n, color, marker in zip(data_lengths, colors, markers):
#     end = start + n

#     ax.scatter(P_all[start:end], T_all[start:end], marker=marker,
#         c=color, s=40)

#     start = end


ax.legend(frameon=False, ncol=2)
ax.grid(True, linestyle=':', color='lightgray', alpha=0.6)

plt.tight_layout()
# plt.savefig('./phase_diagram_h20_nh3.pdf')
plt.show()


########################################
#            Data storage
########################################

# Compute the total Gibbs energy of the liquid mixture
# Use the formulation:
#   Gtot = (1-x)G_H2O + xG_NH3 + G_mix + G_ex. ; x = x_NH3
# Remind: G_mix = RT[(1-x) ln(1-x)+ xln(x)]
# The pure components H2O and NH3 can be easily computed with existing EOS
# The excess Gibbs energy has been computed here using the Redlich-Kister model

P_export = np.linspace(1, 900, 200)
T_export = np.linspace(223, 273, 150)

export_liquid_gibbs(L_fit=L_fit, compositions=compositions,
    P_grid=P_export, T_grid=T_export, C0=0.0)