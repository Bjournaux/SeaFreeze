"""
Build the SeaFreeze 'NH3aq_RK' Gibbs spline from the Redlich-Kister model.

    python SourceData/NH3aq/build_NH3aq_spline.py            (from the repository root)

Writes Python/seafreeze/splines/NH3_aq_RK_2026/NH3_aq_RK_2026.mat and prints the
accuracy of the spline against the analytic model at points between the sites.
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / 'Python'))
sys.path.insert(0, str(HERE))

from seafreeze import splinebuilder as sb          # noqa: E402
from mlbspline.eval import evalMultivarSpline      # noqa: E402
import nh3aq_model as model                        # noqa: E402

OUT = REPO / 'Python' / 'seafreeze' / 'splines' / 'NH3_aq_RK_2026' / 'NH3_aq_RK_2026.mat'

# Domain of the spline
P_MAX = 1000.0      # MPa
T_MIN, T_MAX = 200.0, 310.0
M_MAX = 18.0        # mol/kg (x_NH3 = 0.245); the fit reaches 17.5 mol/kg (x = 0.24)

# Sites: denser at low pressure (fast change of compressibility) and near m = 0
P_SITES = np.unique(np.concatenate([np.linspace(0.0, 100.0, 11), np.linspace(100.0, P_MAX, 37)]))
T_SITES = np.linspace(T_MIN, T_MAX, 67)
M_SITES = sb.dilute_sites(M_MAX, 46, power=3.0)


def main():
    print(f"sites: P {P_SITES.size}, T {T_SITES.size}, m {M_SITES.size}")
    sp = sb.build_gibbs_spline(model.G_solution, P_SITES, T_SITES, M_SITES)

    # ---- accuracy between the sites -------------------------------------------------
    mid = lambda s: 0.5 * (s[:-1] + s[1:])
    Pt, Tt = mid(P_SITES)[::3], mid(T_SITES)[::3]
    mt = np.concatenate([[1e-5, 1e-3, 0.02], mid(M_SITES)[8::4]])
    x = np.empty(3, dtype=object); x[0], x[1], x[2] = Pt, Tt, mt
    G0 = evalMultivarSpline(sp, x, [0, 0, 0]); G1 = evalMultivarSpline(sp, x, [0, 0, 1])
    f = 1.0 + model.M_NH3 * mt
    muw = model.M_H2O * (G0 - f * mt * G1)
    err_muw = np.abs(muw - model.muw_exact(Pt, Tt, mt))
    err_G = np.abs(G0 - model.G_solution(Pt[:, None, None], Tt[None, :, None], mt[None, None, :]))
    print(f"max |muw(spline) - muw(model)| = {err_muw.max():.4f} J/mol "
          f"(~{err_muw.max() / 20:.4f} K on a melting curve)")
    i = np.unravel_index(err_muw.argmax(), err_muw.shape)
    print(f"   at P = {Pt[i[0]]:.1f} MPa, T = {Tt[i[1]]:.1f} K, m = {mt[i[2]]:.3g} mol/kg")
    print(f"max |G(spline) - G(model)|     = {err_G.max():.4f} J/kg")
    if err_muw.max() > 1.0:
        raise RuntimeError("spline does not reproduce the model to 1 J/mol; refine the sites")

    rk = model.RK
    sb.save_gibbs_spline(
        sp, str(OUT), MW=model.M_NH3, nu=1,
        model='Redlich-Kister excess Gibbs energy fitted to ice liquidus data; '
              'L_k = a_k + b_k (T_ref/T - 1) + v_k (P - P_ref)/1000',
        calibrated='muw, aw, ice melting curves (x_NH3 <= 0.24, 9-971 MPa, 210-270 K)',
        not_calibrated='rho, Cp, Cv, vel, Kt, Ks, Kp, alpha, mus, Vm, Cpm, Va, Cpa, Vex and all '
                       'other properties depending on pure NH3 or on derivatives of the excess term',
        water_reference="SeaFreeze water1 for T >= 240 K; constant-Cp continuation below 240 K",
        nh3_reference='CoolProp HEOS::Ammonia liquid, tabulated 1-900 MPa, 233.4-273 K, continued outside',
        rk_names=' '.join(sorted(rk)), rk_values=np.array([rk[k] for k in sorted(rk)], float),
        T_ref=model.T_REF, P_ref=model.P_REF, P_scale=model.P_SCALE,
        fit_rms_K=float(model.FIT_INFO['rms_K']), fit_n_points=float(model.FIT_INFO['n_points']),
        fit_x_range=np.array(model.FIT_INFO['x_range'], float),
        fit_P_range_MPa=np.array(model.FIT_INFO['P_range_MPa'], float),
        fit_T_range_K=np.array(model.FIT_INFO['T_range_K'], float),
        version='NH3_aq_RK_2026')
    print(f"written {OUT}  ({OUT.stat().st_size / 1e6:.2f} MB), coefs {sp['coefs'].shape}")


if __name__ == '__main__':
    main()
