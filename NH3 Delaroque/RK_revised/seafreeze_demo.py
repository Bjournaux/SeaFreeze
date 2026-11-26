"""
Exercise every public SeaFreeze function on the experimental 'NH3aq_RK' material and
plot density / heat-capacity surfaces versus P, T and composition.

    python "NH3 Delaroque/RK_revised/seafreeze_demo.py"      (from the repository root)

Writes outputs/seafreeze_functions_NH3aq.png and outputs/seafreeze_surfaces_NH3aq.png.
Density and Cp are NOT calibrated for NH3aq_RK: the surfaces show what the material
returns, not a validated equation of state.
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "Python"))
import seafreeze as sf                      # noqa: E402
import seafreeze.seafreeze as sfm           # noqa: E402

OUT = HERE / "outputs"
M_NH3, NW = sfm.mNH3_kgmol, 1.0 / sfm.mH2O_kgmol
wt2m = lambda w: (np.asarray(w, float) / (100.0 - np.asarray(w, float))) / M_NH3
m2wt = lambda m: 100.0 * np.asarray(m, float) * M_NH3 / (1.0 + np.asarray(m, float) * M_NH3)
ICE_NAME = {0: "liquid", 1: "Ih", 2: "II", 3: "III", 5: "V", 6: "VI"}
ICE_COL = {"Ih": "#0072B2", "II": "#D55E00", "III": "#009E73", "V": "#CC79A7", "VI": "#E69F00"}

results = []


def check(function, what, ok, detail):
    results.append((function, what, "PASS" if ok else "FAIL", detail))


def grid(P, T, m):
    return np.array([np.atleast_1d(np.asarray(P, float)), np.atleast_1d(np.asarray(T, float)),
                     np.atleast_1d(np.asarray(m, float))], dtype=object)


def exp_warnings(fn):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        out = fn()
    return out, [x for x in w if issubclass(x.category, sf.ExperimentalMaterialWarning)]


# ---------------------------------------------------------------------------
# 1. Registry
# ---------------------------------------------------------------------------
check("is_solution", "NH3aq_RK / NaClaq are solutions, water1 / Ih are not",
      sf.is_solution("NH3aq_RK") and sf.is_solution("NaClaq") and not sf.is_solution("water1") and not sf.is_solution("Ih"),
      "True, True, False, False")
check("solution_phases", "materials per solute",
      sf.solution_phases("NH3") == ["NH3aq_RK"] and len(sf.solution_phases("NaCl")) == 4,
      f"NH3: {sf.solution_phases('NH3')}; NaCl: {len(sf.solution_phases('NaCl'))} materials")
check("solute_phase", "solute name -> default material",
      sf.solute_phase("nh3") == "NH3aq_RK" and sf.solute_phase("NaCl") == "NaClaq", "NH3 -> NH3aq_RK, NaCl -> NaClaq")
rng = sf.phase_range("NH3aq_RK")
check("phase_range", "domain of the NH3aq_RK spline", rng.P == (0.0, 1000.0) and rng.T == (200.0, 310.0) and rng.m == (0.0, 18.0),
      f"P {rng.P} MPa, T {rng.T} K, m {rng.m} mol/kg")
check("uncalibrated_properties", "names what is not calibrated",
      sf.uncalibrated_properties("NH3aq_RK", ("rho", "muw", "Cp")) == ["Cp", "rho"], "['Cp', 'rho'] flagged, muw not")
sf.register_solution("DEMOaq", "DEMO", M_NH3, 1, 1e-6, spline=sfm._SPLINE_MAP["NH3aq_RK"], experimental=True, calibrated=("muw",))
ok = sf.is_solution("DEMOaq") and sf.solute_phase("DEMO") == "DEMOaq"
a = sf.getProp(grid(100, 260, 5), "DEMOaq", sfm.defpath, "muw").muw
b = sf.getProp(grid(100, 260, 5), "NH3aq_RK", sfm.defpath, "muw").muw
sf.unregister_solution("DEMOaq")
check("register_solution / unregister_solution", "a material registered at run time works, then is removed",
      ok and np.array_equal(a, b) and "DEMOaq" not in sf.phases, "muw identical to NH3aq_RK; registry restored")

# ---------------------------------------------------------------------------
# 2. getProp
# ---------------------------------------------------------------------------
PTm = grid([10, 300, 800], [230, 260, 290], [0.0, 2.0, 8.0, 15.0])
allp, w_all = exp_warnings(lambda: sf.getProp(PTm, "NH3aq_RK"))
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    nacl = sf.getProp(grid([10, 300, 800], [250, 280, 300], [0.5, 2.0, 5.0]), "NaClaq")
check("getProp (all properties, grid)", "same property set as NaClaq, all finite, shape (3,3,4)",
      set(vars(allp)) == set(vars(nacl)) and all(np.all(np.isfinite(np.asarray(v, float))) for v in vars(allp).values())
      and allp.rho.shape == (3, 3, 4), f"{len(vars(allp))} properties")
check("getProp (experimental warning)", "warns once when uncalibrated properties are requested", len(w_all) == 1,
      str(w_all[0].message)[:70] + "...")
cal, w_cal = exp_warnings(lambda: sf.getProp(PTm, "NH3aq_RK", sfm.defpath, "muw", "aw", "phi"))
check("getProp (selective)", "returns exactly what was asked, no warning for calibrated properties",
      set(vars(cal)) == {"muw", "aw", "phi"} and not w_cal, "muw, aw, phi; 0 warnings")
sc = np.empty((2,), dtype=object); sc[0], sc[1] = (300.0, 260.0, 8.0), (800.0, 290.0, 15.0)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    s = sf.getProp(sc, "NH3aq_RK", sfm.defpath, "rho", "muw")
check("getProp (scatter)", "scatter points equal the matching grid cells",
      np.allclose(s.rho, [allp.rho[1, 1, 2], allp.rho[2, 2, 3]], rtol=1e-12)
      and np.allclose(s.muw, [allp.muw[1, 1, 2], allp.muw[2, 2, 3]], rtol=1e-12), "rho and muw identical")
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    w1 = sf.getProp(np.array([PTm[0], np.array([260.0, 290.0])], dtype=object), "water1")
check("getProp (pure-water limit)", "m = 0 reproduces water1 (rho, Cp, G) for T >= 240 K",
      np.allclose(allp.rho[:, 1:, 0], w1.rho, rtol=1e-5) and np.allclose(allp.Cp[:, 1:, 0], w1.Cp, rtol=1e-4)
      and np.allclose(allp.G[:, 1:, 0], w1.G, atol=1e-2),
      f"max rel. diff rho {np.max(np.abs(allp.rho[:, 1:, 0] / w1.rho - 1)):.1e}, Cp {np.max(np.abs(allp.Cp[:, 1:, 0] / w1.Cp - 1)):.1e}")
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    o = sf.getProp(grid(5000, 260, 5), "NH3aq_RK", sfm.defpath, "rho", "muw")
check("getProp (out of domain)", "NaN outside the spline domain, no crash", np.all(np.isnan(np.asarray(o.rho, float))) and np.all(np.isnan(np.asarray(o.muw, float))),
      "5000 MPa -> NaN")
aw = np.ravel(sf.getProp(grid(300, 260, [0.5, 2, 5, 10, 15]), "NH3aq_RK", sfm.defpath, "aw").aw)
xw = NW / (NW + np.array([0.5, 2, 5, 10, 15]))
check("getProp (water activity)", "aw decreases with NH3 and falls below the mole fraction of water (negative deviation)",
      np.all(np.diff(aw) < 0) and aw[-1] < xw[-1], f"aw = {np.round(aw, 3)}, x_w = {np.round(xw, 3)}")

# ---------------------------------------------------------------------------
# 3. phase_lines
# ---------------------------------------------------------------------------
m10, m20 = float(wt2m(10)), float(wt2m(20))
lines = {}
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    (_, w_pl) = exp_warnings(lambda: sf.phase_lines("Ih", "NH3aq_RK", m=m10))
    for ice in ("Ih", "II", "III", "V", "VI"):
        lines[ice] = {w: sf.phase_lines(ice, "NH3aq_RK", m=float(wt2m(w))) for w in (5, 10, 20)}
    pure = sf.phase_lines("Ih", "water1")
at = lambda r, P: np.interp(P, r.P[np.argsort(r.P)], r.T[np.argsort(r.P)])
T100 = [at(pure, 100.0)] + [at(lines["Ih"][w], 100.0) for w in (5, 10, 20)]
check("phase_lines (ice-NH3aq_RK)", "melting curve for each of the five ices; no experimental warning",
      all(lines[i][10].P.size > 50 for i in lines) and not w_pl, "Ih, II, III, V, VI")
check("phase_lines (depression)", "ice Ih melting temperature at 100 MPa falls with NH3 content",
      np.all(np.diff(T100) < 0), "0/5/10/20 wt%: " + ", ".join(f"{t:.1f}" for t in T100) + " K")
both = sf.phase_lines("Ih", "NH3aq_RK", m=[m10, m20])
check("phase_lines (vector m, order)", "list of results for several molalities; argument order irrelevant",
      [round(r.m, 3) for r in both] == [round(m10, 3), round(m20, 3)]
      and np.array_equal(sf.phase_lines("NH3aq_RK", "Ih", m=m10).T, lines["Ih"][10].T), "2 results; symmetric")
try:
    sf.phase_lines("Ih", "NH3aq_RK"); ok = False
except ValueError as e:
    ok = "NH3aq_RK" in str(e)
try:
    sf.phase_lines("Ih", "NH3aq_RK", m=50.0); ok2 = False
except ValueError:
    ok2 = True
check("phase_lines (validation)", "missing or out-of-range molality is rejected with a clear error", ok and ok2,
      "m missing -> error naming NH3aq_RK; m = 50 -> error")

# ---------------------------------------------------------------------------
# 4. whichphase
# ---------------------------------------------------------------------------
Pw, Tw = np.arange(0.1, 1000.0, 5.0), np.arange(205.0, 300.0, 0.5)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    wp = {w: np.squeeze(sf.whichphase(grid(Pw, Tw, wt2m(w)), "NH3aq_RK")) for w in (10, 20)}
    wp_name = np.squeeze(sf.whichphase(grid(Pw, Tw, wt2m(10)), "NH3"))


def sequence(arr):
    """Ice just below the liquidus, along pressure."""
    seq = []
    for i in range(arr.shape[0]):
        liq = np.where(arr[i] == 0)[0]
        if liq.size and liq[0] > 0:
            name = ICE_NAME[int(arr[i, liq[0] - 1])]
            if not seq or seq[-1] != name:
                seq.append(name)
    return seq


seq10, seq20 = sequence(wp[10]), sequence(wp[20])
check("whichphase (phase map)", "ice sequence along the liquidus; ice II on the liquidus, no ice III at 10-20 wt%",
      seq10 == ["Ih", "II", "V", "VI"] and seq20 == ["Ih", "II", "V", "VI"], f"10 wt%: {' -> '.join(seq10)}; 20 wt%: {' -> '.join(seq20)}")
check("whichphase (solute name)", "solute='NH3' gives the same map as 'NH3aq_RK'", np.array_equal(wp_name, wp[10], equal_nan=True), "identical")
i100 = np.argmin(np.abs(Pw - 100.1))
T_wp = Tw[np.where(wp[10][i100] == 0)[0][0]]
check("whichphase vs phase_lines", "liquid boundary of the map agrees with the Ih melting line at 100 MPa",
      abs(T_wp - at(lines["Ih"][10], 100.1)) < 0.6, f"map {T_wp:.1f} K, line {at(lines['Ih'][10], 100.1):.2f} K (map step 0.5 K)")
pts = np.empty((2,), dtype=object); pts[0], pts[1] = (100.0, 240.0, m10), (100.0, 270.0, m10)
check("whichphase (scatter)", "ice Ih below, liquid above the liquidus", list(sf.whichphase(pts, "NH3aq_RK")) == [1, 0], "[1, 0]")

# ---------------------------------------------------------------------------
# 5. rho2P
# ---------------------------------------------------------------------------
P_ref = np.array([20.0, 150.0, 400.0, 700.0, 950.0])
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    rho_ref = np.array([sf.getProp(grid(p, 260, 5), "NH3aq_RK", sfm.defpath, "rho").rho.ravel()[0] for p in P_ref])
P_back, w_r = exp_warnings(lambda: sf.rho2P(rho_ref, 260.0, "NH3aq_RK", m=5.0))
check("rho2P", "density -> pressure round trip, and it warns that density is not calibrated",
      np.allclose(P_back, P_ref, atol=0.05) and len(w_r) == 1, f"max error {np.max(np.abs(P_back - P_ref)):.1e} MPa; 1 warning")
try:
    sf.rho2P(rho_ref, 260.0, "NH3aq_RK"); ok = False
except ValueError:
    ok = True
check("rho2P (validation)", "molality is required for a solution", ok, "error raised")

# ---------------------------------------------------------------------------
# 6. wpd + figure of the functions
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(14, 10.5))
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    n0 = len(sf.wpd(show_meta=False).axes[0].get_lines()); plt.close(plt.gcf())   # same options, no solute
    sf.wpd(ax=axes[0, 0], solute="NH3", m=[float(wt2m(w)) for w in (5, 10, 20)], show_meta=False)
ax = axes[0, 0]
check("wpd", "water phase diagram with NH3aq_RK melting curves overlaid", len(ax.get_lines()) > n0,
      f"{len(ax.get_lines()) - n0} NH3aq_RK curves added")
ax.set(xlim=(0, 1000), ylim=(200, 300), title="wpd(solute='NH3', m=5, 10, 20 wt%)")
if ax.get_legend():
    ax.get_legend().remove()

ax = axes[0, 1]
levels = [-0.5, 0.5, 1.5, 2.5, 3.5, 5.5, 6.5]
cols = ["#f2f2f2", ICE_COL["Ih"], ICE_COL["II"], ICE_COL["III"], ICE_COL["V"], ICE_COL["VI"]]
ax.contourf(Pw, Tw, wp[10].T, levels=levels, colors=cols, alpha=0.55)
for ice in lines:
    r = lines[ice][10]
    o = np.argsort(r.P)
    ax.plot(r.P[o], r.T[o], color=ICE_COL[ice], lw=1.6, label=f"phase_lines({ice}, NH3aq_RK)")
for name, (x, y) in {"liquid": (500, 285), "Ih": (70, 215), "II": (330, 215), "V": (540, 232), "VI": (850, 240)}.items():
    ax.text(x, y, name, ha="center", fontsize=10, weight="bold")
ax.set(xlim=(0, 1000), ylim=(205, 300), xlabel="Pressure (MPa)", ylabel="Temperature (K)",
       title="whichphase map at 10 wt% NH3 (colours) with phase_lines (curves)")
ax.legend(fontsize=7, loc="upper left", frameon=False)

ax = axes[1, 0]
mm = np.linspace(0, 18, 73)
for (P, T), c in zip(((0.1, 260), (300, 250), (800, 270)), ("#0072B2", "#D55E00", "#009E73")):
    a_ = np.ravel(sf.getProp(grid(P, T, mm), "NH3aq_RK", sfm.defpath, "aw").aw)
    ax.plot(m2wt(mm), a_, color=c, lw=1.8, label=f"aw, {P:g} MPa, {T} K")
ax.plot(m2wt(mm), NW / (NW + mm), "k--", lw=1, label="ideal solution (mole fraction of water)")
ax.set(xlabel="NH3 (wt%)", ylabel="Water activity", title="getProp(..., 'aw'): calibrated quantity")
ax.legend(fontsize=8, frameon=False); ax.grid(True, ls=":", c="lightgray")

ax = axes[1, 1]
ax.plot(P_ref, P_back - P_ref, "o-", color="#0072B2")
ax.axhline(0, c="k", lw=0.8)
ax.set(xlabel="Pressure used to compute the density (MPa)", ylabel="rho2P result - pressure (MPa)",
       title="rho2P round trip at 260 K, 5 mol/kg (density not calibrated)", ylim=(-0.05, 0.05))
ax.grid(True, ls=":", c="lightgray")
fig.suptitle("SeaFreeze functions on the experimental material 'NH3aq_RK'", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(OUT / "seafreeze_functions_NH3aq.png", dpi=140)
plt.close(fig)

# ---------------------------------------------------------------------------
# 7. Surfaces of density and Cp versus P, T, composition
# ---------------------------------------------------------------------------
Pg, Tg, wg = np.linspace(1.0, 1000.0, 41), np.linspace(210.0, 300.0, 37), np.linspace(0.0, 23.0, 36)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    full = sf.getProp(grid(Pg, Tg, wt2m(wg)), "NH3aq_RK", sfm.defpath, "rho", "Cp")
    liq = np.stack([np.squeeze(sf.whichphase(grid(Pg, Tg, m_), "NH3aq_RK")) == 0 for m_ in wt2m(wg)], axis=2)
check("getProp (large 3-D grid)", "41 x 37 x 36 grid evaluates, finite everywhere",
      full.rho.shape == (41, 37, 36) and np.all(np.isfinite(full.rho)) and np.all(np.isfinite(full.Cp)),
      f"rho {full.rho.min():.0f}-{full.rho.max():.0f} kg/m3, Cp {full.Cp.min():.0f}-{full.Cp.max():.0f} J/kg/K")

iw, iT, iP = np.argmin(np.abs(wg - 10.0)), np.argmin(np.abs(Tg - 260.0)), np.argmin(np.abs(Pg - 300.0))
slices = [
    (Pg, Tg, lambda a: a[:, :, iw], "Pressure (MPa)", "Temperature (K)", f"{wg[iw]:.0f} wt% NH3"),
    (Pg, wg, lambda a: a[:, iT, :], "Pressure (MPa)", "NH3 (wt%)", f"{Tg[iT]:.0f} K"),
    (Tg, wg, lambda a: a[iP, :, :], "Temperature (K)", "NH3 (wt%)", f"{Pg[iP]:.0f} MPa"),
]
fig = plt.figure(figsize=(17, 10.5))
for r, (prop, label, cmap) in enumerate((("rho", "Density (kg/m$^3$)", "viridis"), ("Cp", "Cp (J/kg/K)", "magma"))):
    for c, (xa, ya, cut, xl, yl, ttl) in enumerate(slices):
        ax = fig.add_subplot(2, 3, r * 3 + c + 1, projection="3d")
        Z = cut(getattr(full, prop)); L = cut(liq)
        X, Y = np.meshgrid(xa, ya, indexing="ij")
        norm = plt.Normalize(np.nanmin(Z), np.nanmax(Z))
        colors = plt.get_cmap(cmap)(norm(Z))
        colors[~L] = (0.75, 0.75, 0.75, 1.0)             # liquid metastable: ice is the stable phase
        ax.plot_surface(X, Y, Z, facecolors=colors, rstride=1, cstride=1, linewidth=0, antialiased=True, shade=False)
        ax.set_xlabel(xl, fontsize=8); ax.set_ylabel(yl, fontsize=8); ax.set_zlabel(label, fontsize=8)
        ax.set_title(f"{label.split(' (')[0]} at {ttl}", fontsize=10)
        ax.tick_params(labelsize=7); ax.view_init(elev=24, azim=-58)
fig.suptitle("getProp('NH3aq_RK'): density and heat capacity versus pressure, temperature and composition\n"
             "NOT calibrated for this material. Colour: value where the liquid is stable; grey: liquid metastable (ice stable, from whichphase)",
             fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(OUT / "seafreeze_surfaces_NH3aq.png", dpi=140)
plt.close(fig)

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
w0 = max(len(r[0]) for r in results); w1 = max(len(r[1]) for r in results)
print(f"\n{'function':<{w0}}  {'check':<{w1}}  result  detail")
print("-" * (w0 + w1 + 30))
for f_, what, res, det in results:
    print(f"{f_:<{w0}}  {what:<{w1}}  {res:<6}  {det}")
n_fail = sum(r[2] == "FAIL" for r in results)
print(f"\n{len(results) - n_fail} of {len(results)} checks passed")
print("figures:", OUT / "seafreeze_functions_NH3aq.png", "and", OUT / "seafreeze_surfaces_NH3aq.png")
sys.exit(1 if n_fail else 0)
