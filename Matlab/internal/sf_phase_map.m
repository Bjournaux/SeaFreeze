function pm = sf_phase_map(P, T, varargin)
% SF_PHASE_MAP  Stable phase of H2O on a (P,T) grid by Gibbs-energy minimisation.
%
%   pm = sf_phase_map(P, T)
%   pm = sf_phase_map(P, T, 'fluid', 'water_Brown2026', 'ices', {'Ih','II','III','V','VI'}, ...
%                     'dilute_extension', true, 'sanity', true, 'melt_mask', true)
%
% The fluid is a Helmholtz material (vapour, liquid, supercritical); the
% ices are SeaFreeze Gibbs splines.  All share the IAPWS-95 reference state,
% so G per kg is compared directly.
%   dilute_extension  below the fluid's lowest T (230 K), use its ideal-gas
%                     part at P < 1e-4 MPa (vapour field / sublimation)
%   sanity            an ice competes only where its spline is physical
%                     (rho > 0, Kt > 0, 0 < Cp < 2 x 9R/M)
%   melt_mask         the fluid is ignored more than 40 K below the melting
%                     curve of the stable solid (sf_melt_T_dq2026) above the
%                     triple-point pressure — the surface's validity mask
%
% Output struct: P (nP x 1), T (1 x nT), names {fluid, ices...},
%   G, rho (nphase x nP x nT), stable (nP x nT, index into names, 0 = none),
%   rho_stable (nP x nT).
p = inputParser;
addParameter(p, 'fluid', 'water_Brown2026');
addParameter(p, 'ices', {'Ih','II','III','V','VI'});
addParameter(p, 'dilute_extension', true);
addParameter(p, 'sanity', true);
addParameter(p, 'melt_mask', true);
parse(p, varargin{:});
fluid = char(p.Results.fluid); ices = cellstr(p.Results.ices);
P = P(:); T = T(:).'; nP = numel(P); nT = numel(T);
names = [{fluid}, ices]; n = numel(names);
G = NaN(n, nP, nT); rho = G;

o = SF_getprop({P, T}, fluid, {'G', 'rho'});
Gf = o.G; rf = o.rho;
if p.Results.melt_mask
    Tm = sf_melt_T_dq2026(P);
    bad = bsxfun(@and, P > 611.657e-6, bsxfun(@lt, T, Tm - 40));
    Gf(bad) = NaN;
end
if p.Results.dilute_extension
    sp = sf_load_spline(fluid);
    Tmin = sp.Tc * exp(sp.knots{2}(1));
    [Pm, Tm2] = ndgrid(P, T);
    lo = Tm2 < Tmin & Pm < 1e-4;
    if any(lo(:))
        [Gf(lo), rf(lo)] = sf_ideal_gas(sp, Pm(lo), Tm2(lo));
    end
end
G(1, :, :) = Gf; rho(1, :, :) = rf;

cpmax = 2 * 9 * 8.314462618 / 0.018015268;
for k = 1:numel(ices)
    o = SF_getprop({P, T}, ices{k}, {'G', 'rho', 'Cp', 'Kt'});
    g = o.G; g(g == 0) = NaN;                       % out-of-range sentinel
    if p.Results.sanity
        g(~(o.rho > 0 & o.Kt > 0 & o.Cp > 0 & o.Cp < cpmax)) = NaN;
    end
    G(k+1, :, :) = g; rho(k+1, :, :) = o.rho;
end

Gs = G; Gs(isnan(Gs)) = Inf;
[gmin, stable] = min(Gs, [], 1);
stable = reshape(stable, nP, nT); gmin = reshape(gmin, nP, nT);
stable(~isfinite(gmin)) = 0;
rs = NaN(nP, nT);
for k = 1:n
    m = stable == k; rk = reshape(rho(k, :, :), nP, nT); rs(m) = rk(m);
end
pm = struct('P', P, 'T', T, 'names', {names}, 'G', G, 'rho', rho, 'stable', stable, 'rho_stable', rs);
end
