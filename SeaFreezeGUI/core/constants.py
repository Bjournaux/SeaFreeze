"""Material definitions and property metadata for the SeaFreeze GUI."""

# ── Material groups (display order) ──────────────────────────────────────────
MATERIALS_ICE = ["Ih", "II", "III", "V", "VI", "VII_X_French"]
MATERIALS_WATER = ["water1", "water2", "water_IAPWS95"]
MATERIALS_NACL = ["NaClaq", "NaClaq_LP", "NaClaq_HP", "NaClaq_5GPa_2024"]
ALL_MATERIALS = MATERIALS_ICE + MATERIALS_WATER + MATERIALS_NACL

# Human-readable labels — full citation style (used in dropdowns / selectors)
MATERIAL_LABELS = {
    "Ih":               "Ice Ih (Feistel & Wagner, 2006)",
    "II":               "Ice II (Journaux et al., 2020)",
    "III":              "Ice III (Journaux et al., 2020)",
    "V":                "Ice V (Journaux et al., 2020)",
    "VI":               "Ice VI (Journaux et al., 2020)",
    "VII_X_French":     "Ice VII/X (French & Redmer, 2015)",
    "water1":           "Water (Bollengier et al., 2019)",
    "water2":           "Water (Brown, 2018)",
    "water_IAPWS95":    "Water IAPWS95 (Wagner & Pruss, 2002)",
    "NaClaq":           "NaCl(aq) (Brown et al., under review)",
    "NaClaq_LP":        "NaCl(aq) LP (Brown et al., under review)",
    "NaClaq_HP":        "NaCl(aq) HP (Brown et al., under review)",
    "NaClaq_5GPa_2024": "NaCl(aq) 5 GPa (Brown, 2024)",
}

# Short labels — for plot legends and phase diagram traces
MATERIAL_SHORT_LABELS = {
    "Ih":               "Ice Ih",
    "II":               "Ice II",
    "III":              "Ice III",
    "V":                "Ice V",
    "VI":               "Ice VI",
    "VII_X_French":     "Ice VII/X",
    "water1":           "Water",
    "water2":           "Water (Brown)",
    "water_IAPWS95":    "Water IAPWS95",
    "NaClaq":           "NaCl(aq)",
    "NaClaq_LP":        "NaCl(aq) LP",
    "NaClaq_HP":        "NaCl(aq) HP",
    "NaClaq_5GPa_2024": "NaCl(aq) 5GPa",
}

def is_nacl(material: str) -> bool:
    return material in MATERIALS_NACL

def is_solid(material: str) -> bool:
    return material in MATERIALS_ICE

# ── Property metadata: (display name, unit) ──────────────────────────────────
# Ordered as they should appear in tables / selectors.
PROPS_ALL = {
    "G":     ("Gibbs energy",                "J/kg"),
    "S":     ("Entropy",                     "J/K/kg"),
    "U":     ("Internal energy",             "J/kg"),
    "H":     ("Enthalpy",                    "J/kg"),
    "A":     ("Helmholtz energy",            "J/kg"),
    "rho":   ("Density",                     "kg/m³"),
    "Cp":    ("Heat capacity (const P)",     "J/kg/K"),
    "Cv":    ("Heat capacity (const V)",     "J/kg/K"),
    "Kt":    ("Isothermal bulk modulus",     "MPa"),
    "Kp":    ("Pressure derivative of K",    "—"),
    "Ks":    ("Isentropic bulk modulus",     "MPa"),
    "alpha": ("Thermal expansivity",         "1/K"),
    "vel":   ("Sound velocity",              "m/s"),
    "Js":    ("Joule–Thomson coeff.",   "K/MPa"),
    "gamma_Gruneisen": ("Grüneisen parameter", "—"),
}

PROPS_SOLID = {
    "shear": ("Shear modulus",    "MPa"),
    "Vp":    ("P-wave velocity",  "m/s"),
    "Vs":    ("S-wave velocity",  "m/s"),
}

PROPS_NACL = {
    "mus":  ("Solute chemical potential",  "J/mol"),
    "muw":  ("Solvent chemical potential", "J/mol"),
    "Vm":   ("Partial molar vol. solute",  "cm³/mol"),
    "Vw":   ("Partial molar vol. solvent", "cm³/mol"),
    "Va":   ("Apparent molar volume",      "cm³/mol"),
    "Cpa":  ("Apparent molar Cp",          "J/mol/K"),
    "Cpm":  ("Partial molar Cp",           "J/mol/K"),
    "Vex":  ("Excess volume",              "cm³/mol"),
    "phi":  ("Osmotic coefficient",        "—"),
    "aw":   ("Water activity",             "—"),
}


PROP_CATEGORIES = [
    ("Thermodynamic potentials", ["G", "S", "U", "H", "A"]),
    ("Elastic & volumetric",     ["rho", "alpha", "Kt", "Kp", "Ks",
                                  "vel", "shear", "Vp", "Vs"]),
    ("Thermal",                  ["Cp", "Cv", "Js", "gamma_Gruneisen"]),
    ("Mixing (NaCl)",            ["mus", "muw", "Vm", "Vw", "Va",
                                  "Cpa", "Cpm", "Vex", "phi", "aw"]),
]


def available_properties(material: str) -> dict:
    """Return {symbol: (name, unit)} for the given material."""
    props = dict(PROPS_ALL)
    if is_solid(material):
        props.update(PROPS_SOLID)
    if is_nacl(material):
        props.update(PROPS_NACL)
    return props


def categorized_properties(material: str) -> list:
    """Return [(category_name, [(symbol, display_name, unit), ...])] for material.

    Only includes categories that have at least one property available.
    """
    avail = available_properties(material)
    categories = []
    for cat_name, symbols in PROP_CATEGORIES:
        items = [(s, avail[s][0], avail[s][1]) for s in symbols if s in avail]
        if items:
            categories.append((cat_name, items))
    return categories
