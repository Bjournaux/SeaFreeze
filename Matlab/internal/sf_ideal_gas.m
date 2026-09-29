function [G, rho] = sf_ideal_gas(sp, P, T)
% SF_IDEAL_GAS  Gibbs energy (J/kg) and density (kg/m^3) of a psi surface's
% ideal-gas part alone (Z = 1) at P (MPa), T (K).
%
% This is the dilute-vapour limit of the EOS.  It needs no spline and so
% stays defined below the surface's lowest temperature knot; used by
% SF_coexistence ('dilute_extension') and sf_phase_map.
R = double(sp.R); rhoc = double(sp.rhoc); Tc = double(sp.Tc);
P = P(:); T = T(:);
rho = P * 1e6 ./ (R * T);
d = rho / rhoc; tau = Tc ./ T;
n0 = double(sp.phi0_n0(:)'); g0 = double(sp.phi0_g0(:)');
phi0 = log(d) + n0(1) + n0(2) * tau + n0(3) * log(tau);
for i = 1:numel(g0)
    phi0 = phi0 + n0(3+i) * log(1 - exp(-g0(i) * tau));
end
G = R * T .* (phi0 + 1);
end
