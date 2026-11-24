from pathlib import Path
import numpy as np
import pandas as pd

from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import brentq, least_squares

from CoolProp.CoolProp import PropsSI

# Physical constants
R = 8.314462618       # in J/mol/K, gas constant
M_H2O = 0.01801528    # converstion factor J/kg->J/mol

P_REF = 0.1           # MPa, surface pressure


##############################
# 1.      Data loading
##############################

def read_seafreeze_table(filename):

    with open(filename, "r") as f:
        lines = f.readlines()

    temp_line = next((line for line in lines if "Temperature (K)" in line), None)
    pres_line = next((line for line in lines if "Pressure (MPa)" in line), None)
    data_start = next((i + 1 for i, line in enumerate(lines) if "####" in line),  None)

    if temp_line is None or pres_line is None or data_start is None:
        raise ValueError(f"Invalid SeaFreeze file: {filename}")

    # Temperature grid [K]
    T_start, T_end, nT = map(float, temp_line.split("Temperature (K)")[1].strip().split(","))

    T = np.linspace(T_start, T_end, int(nT))

    # Pressure grid [MPa]
    P_start, P_end, nP = map(float, pres_line.split("Pressure (MPa)")[1].strip().split(","))

    P = np.linspace(P_start, P_end, int(nP))

    # Gibbs energy data
    flat = np.fromstring(",".join(lines[data_start:]), sep=",")

    G = flat.reshape(len(P), len(T))

    return P, T, G


def load_liquidus(filename, x_value):

    data = np.loadtxt(filename, delimiter=",", skiprows=1)

    P = data[:, 0]
    T = data[:, 1]
    x = np.full(P.shape, x_value)

    return P, T, x

# Data from Choukroun & Grasset 2010, fig. 1 and Leliwa 2002
# acquired for solution of several molar fraction x_NH3 values
def load_liquidus_datasets(datasets):

    P_all = []
    T_all = []
    x_all = []

    for filename, x_value in datasets:

        P, T, x = load_liquidus(filename, x_value)

        P_all.append(P)
        T_all.append(T)
        x_all.append(x)

    return (np.concatenate(P_all), np.concatenate(T_all), np.concatenate(x_all))



##############################
# 2.     Thermodynamics
##############################

# Here, fitting data and find A0 
# Table for G_liq_H20 and G_Ih
def initialize_fusion_model(water_file, ice_file):

    P_liq, T_liq, G_liq = read_seafreeze_table(water_file)
    P_ih, T_ih, G_ih = read_seafreeze_table(ice_file)

    # Conversion J/kg -> J/mol
    G_liq *= M_H2O
    G_ih *= M_H2O

    # Interpolate points in the tables of seafereze
    interp_liq = RegularGridInterpolator((P_liq, T_liq),  G_liq)
    interp_ih = RegularGridInterpolator((P_ih, T_ih), G_ih)

    return interp_liq, interp_ih

# Compute ΔG_fus(P,T) = G_Ih_H20 - G_sol_H20
def deltaG_fusion(P, T, interp_liq, interp_ih):

    P, T = np.broadcast_arrays(P, T)

    original_shape = P.shape

    points = np.column_stack([P.ravel(), T.ravel()])

    # Clip independently to each interpolation grid
    points_liq = points.copy()
    points_ih = points.copy()

    for points_i, interp in [(points_liq, interp_liq), (points_ih, interp_ih)]:

        P_grid, T_grid = interp.grid

        points_i[:, 0] = np.clip(points_i[:, 0], P_grid.min(), P_grid.max())
        points_i[:, 1] = np.clip(points_i[:, 1], T_grid.min(), T_grid.max())

    result = interp_ih(points_ih) - interp_liq(points_liq)

    return result.reshape(original_shape)



##############################
# 3.    Activity coeff
##############################

# Activity coefficient from eq. 10.4-6 and 10.4-7 
# => Redlich-Kister model (i.e., 3 free parameters L_i)
def ln_gamma_RK_model(x2, L0, L1, L2):
    x1 = 1.0 - x2
    return x2**2 * (L0 + L1 * (4*x1 - 1) + L2 * (12*x1**2 - 8*x1 + 1))

def RK_coeffs(P, L_fit, C0=0.0, C1=0.0, C2=0.0):

    L0_fit, L1_fit, L2_fit = L_fit

    dP = P - P_REF

    L0 = L0_fit + C0 * dP
    L1 = L1_fit + C1 * dP
    L2 = L2_fit + C2 * dP

    return L0, L1, L2

def ln_gamma_RK(P, x2, L_fit, C0=0.0, C1=0.0, C2=0.0):
    L0, L1, L2 = RK_coeffs(P, L_fit, C0, C1, C2)

    return ln_gamma_RK_model(x2, L0, L1, L2)

##############################
# 4. Redlich-Kister fitting
##############################

# For each experimental points (P,T,x2), compute the residual between what
# the thermodynamical model's predicts and what the ice-liq eq. dictates.
def residuals(params, P, T, x2, interp_liq, interp_ih):
    L0, L1, L2 = params
    x1 = 1.0 - x2

    lngamma = ln_gamma_RK_model(x2, L0, L1, L2)

    # NOTA BENE -- FEW EXPLANATIONS
    # The liquid chemical potential of a non-ideal solution:     μ1_liq = μ1_liqº + RT ln(a1) ; 
    # where º stands for pure/reference and a1 = x1 * γ1.
    # At the equilibrium ice-liquid, the condition becomes:      μ_Ih_H20 = μ_sol_H20, 
    # where the chem. pot. of the solution is written as:        μ_sol_H20 = μ_liqº + RT*ln(a1).
    # This eq. is equivalent to:                                 G_Ih_H20 = G_liqº + RT*ln(a1)
    # OR G_Ih_H20 - G_liqº = RT*ln(a1)
    # As G_Ih_H20 - G_liqº = ΔG_fus, then the eq. becomes ΔG_fus = RT * ln(a1)
    # The residual is equivalent to ΔG_fus - RT ln(a1). 
    model = R * T * (np.log(x1) + lngamma)
    obs = deltaG_fusion(P, T, interp_liq, interp_ih)

    return (model - obs) / (R * T)

def residuals_V(params, P, Tobs, x, interp_liq, interp_V, L_fit):
    # C0 parameters is calibrated obly with values of experimental datapoint of ice V
    # the datapoints of ice III are less reliable...
    C0 = params[0]

    x1 = 1.0 - x

    lngamma = ln_gamma_RK(P, x, L_fit, C0=C0)
    model = R * Tobs * (np.log(x1) + lngamma)
    obs = deltaG_solid_liquid(P, Tobs, interp_liq, interp_V)

    return (model - obs) / (R * Tobs)

def fit_margules(P, T, x, interp_liq, interp_ih):

    guess = [0.0, 0.0, 0.0]

    fit = least_squares(
        residuals,
        guess,
        bounds=([-20, -20, -20], [20, 20, 20]),
        args=(P, T, x, interp_liq, interp_ih)
    )

    if not fit.success:
        raise RuntimeError(f"Margules fitting failed: {fit.message}")

    return fit

def deltaG_solid_liquid(P, T, interp_liq, interp_solid):

    P, T = np.broadcast_arrays(
        np.asarray(P, dtype=float),
        np.asarray(T, dtype=float)
    )

    points = np.column_stack([P.ravel(), T.ravel()])

    def interpolate(interp):

        pts = points.copy()
        P_grid, T_grid = interp.grid

        pts[:, 0] = np.clip(pts[:, 0], P_grid.min(), P_grid.max())
        pts[:, 1] = np.clip(pts[:, 1], T_grid.min(), T_grid.max())

        return interp(pts)

    result = interpolate(interp_solid) - interpolate(interp_liq)

    return result.reshape(P.shape)


##############################
# 5. Liquidus calculations
##############################

def liquidus_temperature(P, x2, L0, L1, L2, interp_liq, interp_ih, T_bounds=(180, 280)):

    x1 = 1.0 - x2

    lngamma = ln_gamma_RK_model(x2, L0, L1, L2)

    rhs_factor = R * (np.log(x1) + lngamma)

    def f(T):
        lhs = deltaG_fusion(P, T, interp_liq, interp_ih)
        return lhs - T * rhs_factor

    return brentq(f, *T_bounds)



def find_liquidus_temperature(P, x2, interp_liq, interp_solid, L_fit,
    T_bounds=(180, 300), n_grid=120, C0=0.0, C1=0.0, C2=0.0,):

    """
    Find the liquidus temperature for a given pressure and NH3 mole fraction.
    """

    if not 0.0 <= x2 < 1.0:
        raise ValueError("x2 must be in the interval [0, 1).")

    if x2 == 0.0:
        rhs_factor = 0.0
    else:
        if L_fit is None:
            raise ValueError("L_fit must be provided for x2 > 0.")

        lngamma = ln_gamma_RK(P, x2, L_fit, C0=C0, C1=C1, C2=C2)

        rhs_factor = R * (np.log(1.0 - x2) + lngamma)

    T_grid = np.linspace(*T_bounds, n_grid)

    deltaG = deltaG_solid_liquid(np.full_like(T_grid, P), T_grid,
        interp_liq, interp_solid)

    residual = deltaG - T_grid * rhs_factor

    # Check for an exact root on the grid
    exact_roots = np.where(np.isclose(residual, 0.0, atol=1e-10))[0]

    if len(exact_roots):
        return T_grid[exact_roots[0]]

    # Find sign changes
    crossings = np.where(np.signbit(residual[:-1]) != np.signbit(residual[1:]))[0]

    if len(crossings) == 0:
        return np.nan

    i = crossings[0]

    def f(T):
        return (deltaG_solid_liquid(P, T, interp_liq, interp_solid) - T * rhs_factor)

    return brentq(f, T_grid[i], T_grid[i + 1])

def liquidus_global(P, x, ice_models, L_fit):

    solutions = []

    for phase, model in ice_models.items():

        T_liq = find_liquidus_temperature(P=P, x2=x, interp_liq=model["liquid"], interp_solid=model["solid"],
            L_fit=L_fit, C0=model["C0"], C1=model["C1"], C2=model["C2"], T_bounds=(180, 300),)

        if np.isfinite(T_liq):
            solutions.append((phase, T_liq))

    if not solutions:
        return np.nan, None

    return max(solutions, key=lambda item: item[1])[::-1]


def calculate_liquidus_curves(P_range, x_values, interp_liq,
    interp_solid, L_fit, T_bounds=(180, 300),  C0=0.0, C1=0.0, C2=0.0):

    curves = {}

    for x in x_values:

        curves[x] = np.array([find_liquidus_temperature(P, x, interp_liq, interp_solid, L_fit=L_fit, 
                T_bounds=T_bounds, C0=C0, C1=C1, C2=C2) for P in P_range])

    return curves


##############################
# 6. Phase boundary
##############################

def load_phase_boundary(filename):
    """
    Load a pure-water phase boundary.
    File columns: Temperature [K], Pressure [MPa].
    """
    data = np.loadtxt(filename, skiprows=3)

    T = data[:, 0]
    P = data[:, 1]

    return T, P


def load_triple_points(filename):
    """
    Load triple points from a file containing successive
    blocks of phase label, temperature and pressure.
    """
    triples = []

    with open(filename) as f:
        lines = [line.strip() for line in f if line.strip()]

    i = 3  # Skip header

    while i < len(lines):

        if i + 2 >= len(lines):
            raise ValueError(f"Incomplete triple-point entry in {filename}")

        triples.append({"name": lines[i], "T": float(lines[i + 1]), "P": float(lines[i + 2]),})

        i += 3

    return triples


def crop_transition(T, P, triple_points):

    if len(triple_points) == 0:
        return T, P

    T_limits = sorted(tp["T"] for tp in triple_points)

    if len(T_limits) == 1:
        mask = T <= T_limits[0]
    else:
        mask = (T >= T_limits[0]) & (T <= T_limits[-1])

    return T[mask], P[mask]


def interp_transition(T, Tcurve, Pcurve):

    Tcurve = np.asarray(Tcurve)
    Pcurve = np.asarray(Pcurve)

    valid = np.isfinite(Tcurve) & np.isfinite(Pcurve)

    Tcurve = Tcurve[valid]
    Pcurve = Pcurve[valid]

    order = np.argsort(Tcurve)

    return np.interp(T, Tcurve[order], Pcurve[order], left=np.nan, right=np.nan)

def stable_phase(T, P, phase_boundaries):

    def boundary_pressure(name, T):
        boundary = phase_boundaries[name]

        return interp_transition(
            T,
            boundary["T"],
            boundary["P"]
        )

    P_Ih_III = boundary_pressure("Ih_III", T)
    P_III_V = boundary_pressure("III_V", T)
    P_II_III = boundary_pressure("II_III", T)
    P_II_V = boundary_pressure("II_V", T)
    P_II_VI = boundary_pressure("II_VI", T)

    # This is only a placeholder for the phase-region logic.
    # The conditions must be checked against the full phase diagram.

    if np.isfinite(P_Ih_III) and P < P_Ih_III:
        return "Ih"

    # Complete the phase-region logic using the topology
    # of the pure-water phase diagram.

    return None


##############################
# 7. Data storage
##############################

# Output directory
OUTPUT_DIR = Path("./outputs")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def export_gibbs_mixing(L_fit, compositions, P_grid, T_grid,
    C0=0.0, output_dir="./outputs",):

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    L0_fit, L1_fit, L2_fit = L_fit

    P_mesh, T_mesh = np.meshgrid(P_grid, T_grid, indexing="ij")

    all_tables = []

    for x2 in compositions:

        x1 = 1.0 - x2

        # Pressure-dependent RK coefficients
        dP = P_mesh - P_REF

        L0 = L0_fit + C0 * dP
        L1 = L1_fit
        L2 = L2_fit

        # Ideal mixing contribution
        G_ideal = R * T_mesh * (x1 * np.log(x1) + x2 * np.log(x2))

        # Excess Gibbs energy
        G_excess = R * T_mesh * x1 * x2 * (L0 + L1 * (x1 - x2) + L2 * (x1 - x2)**2)

        # Total Gibbs energy of mixing
        G_mix = G_ideal + G_excess

        df = pd.DataFrame({"Pressure_MPa": P_mesh.ravel(), "Temperature_K": T_mesh.ravel(), "x_NH3": x2,
            "G_ideal_J_mol": G_ideal.ravel(), "G_excess_J_mol": G_excess.ravel(), "G_mixing_J_mol": G_mix.ravel(),})

        filename = output_dir / f"Gibbs_mixing_x{x2:.4f}.csv"

        df.to_csv(filename, index=False)

        all_tables.append(df)

        print(f"Saved: {filename}")

    # Also save a single file containing all compositions
    df_all = pd.concat(all_tables, ignore_index=True)
    df_all.to_csv(output_dir / "Gibbs_mixing_all_compositions.csv", index=False)

    print("Saved: Gibbs_mixing_all_compositions.csv")



def calculate_liquid_gibbs(P, T, x, L_fit, C0=0.0, C1=0.0, C2=0.0,):
    """
    Calculate molar Gibbs energy of liquid H2O-NH3 mixture.

    Parameters
    ----------
    P : float
        Pressure [MPa].
    T : float
        Temperature [K].
    x : float
        Molar fraction of NH3.
    L_fit : array
        Redlich-Kister parameters [L0, L1, L2].
    C0, C1, C2 : float
        Pressure-dependence coefficients.
        Pressure is expressed in MPa.

    Returns
    -------
    dict
        Gibbs energy contributions [J/mol].
    """

    x2 = x
    x1 = 1.0 - x2

    # Convert pressure to Pa for CoolProp
    P_Pa = P * 1e6

    # Pure-component Gibbs energies
    G_H2O = PropsSI("Gmolar", "P", P_Pa, "T|liquid", T, "HEOS::Water")
    G_NH3 = PropsSI("Gmolar", "P", P_Pa, "T|liquid", T,  "HEOS::Ammonia")

    # Pure-component contribution
    G_pure = x1 * G_H2O + x2 * G_NH3

    # Ideal mixing
    G_mix = R * T * (x1 * np.log(x1) + x2 * np.log(x2))

    # Pressure-dependent RK parameters
    L0_fit, L1_fit, L2_fit = L_fit

    dP = P - P_REF

    L0 = L0_fit + C0 * dP
    L1 = L1_fit + C1 * dP
    L2 = L2_fit + C2 * dP

    # Excess Gibbs energy
    G_ex = R * T * x1 * x2 * (L0 + L1 * (x1 - x2) + L2 * (x1 - x2)**2)

    # Total Gibbs energy
    G_total = G_pure + G_mix + G_ex

    return {"G_H2O_pure": G_H2O, "G_NH3_pure": G_NH3, "G_pure": G_pure,
        "G_mix": G_mix,  "G_ex": G_ex, "G_liquid": G_total,}

def export_liquid_gibbs(L_fit, compositions, P_grid, T_grid,
    C0=0.0, C1=0.0, C2=0.0, output_dir="./outputs",):

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for x in compositions:

        rows = []

        for P in P_grid:
            for T in T_grid:

                try:
                    result = calculate_liquid_gibbs(
                        P, T, x, L_fit,
                        C0=C0,
                        C1=C1,
                        C2=C2,
                    )

                    rows.append({"Pressure_MPa": P, "Temperature_K": T, "x_NH3": x, **result})

                except ValueError:
                    continue

        df = pd.DataFrame(rows)

        filename = output_dir / f"Gliq_x{x:.4f}.csv"

        df.to_csv(filename, index=False)

        print(f"Saved: {filename}")