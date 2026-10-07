function out = SF_coexistence(kind, T, varargin)
% SF_coexistence  Saturation and sublimation curves from a Helmholtz fluid.
% Baptiste Journaux - 2026
%
% A Helmholtz fluid (water_Brown2026) spans liquid and vapour, so SeaFreeze can solve
%   saturation:   G_liquid(P,T) = G_vapour(P,T)
%   sublimation:  G_ice(P,T)    = G_vapour(P,T)
% per temperature by Newton iteration in ln P:
%   d(G_A - G_B)/d ln P = P (1/rho_A - 1/rho_B).
% Ice Gibbs splines and water_Brown2026 share the IAPWS-95 reference state (U = S = 0
% for liquid at the triple point), so G is directly comparable (J/kg).
%
% Usage:
%   out = SF_coexistence('saturation',  T)                 % T < Tc, water_Brown2026
%   out = SF_coexistence('sublimation', T)                 % ice Ih + water_Brown2026 vapour
%   out = SF_coexistence('sublimation', T, 'ice', 'Ih', 'fluid', 'water_Brown2026')
%   out = SF_coexistence('sublimation', T, 'dilute_extension', false)
%
% Output struct (column vectors, one entry per T):
%   out.T      (K)
%   out.P      (MPa)       coexistence pressure; NaN where not found
%   out.rho_A  (kg/m^3)    liquid (saturation) or ice (sublimation)
%   out.rho_B  (kg/m^3)    vapour
%
% 'dilute_extension' (sublimation only, default true): below the fluid's
%   lowest temperature (230 K for water_Brown2026) treat the vapour as the surface's
%   ideal-gas part alone (Z = 1).  At those sublimation pressures (< 10 Pa)
%   the neglected terms change p_sub by ~1e-5 relative (validated against
%   the NIST measurements of Bielska et al. 2013, 175-253 K).  It is an
%   extrapolation of the fluid surface: a warning is issued the first time
%   it is used in a session (id 'SeaFreeze:diluteExtension').  Pass false to
%   get NaN below the surface instead.
%
% Starting values come from the IAPWS-95 auxiliary vapour-pressure equation
% (Wagner & Pruss 2002, eq. 2.5) and IAPWS R14-08 (Wagner et al. 2011); the
% result is the model's own coexistence.
%
% See also: SF_getprop, SF_PhaseLines.

p = inputParser;
addParameter(p, 'ice', 'Ih', @(s) ischar(s) || isstring(s));
addParameter(p, 'fluid', 'water_Brown2026', @(s) ischar(s) || isstring(s));
addParameter(p, 'dilute_extension', true, @(x) islogical(x) || isnumeric(x));
parse(p, varargin{:});
ice = char(p.Results.ice); fluid = sf_material_name(char(p.Results.fluid));   % renamed materials
dilute = logical(p.Results.dilute_extension);

defs = sf_material_defs();
if ~ismember(fluid, defs.helmholtz_phases)
    error('SF_coexistence:badInput', 'fluid must be a Helmholtz material (%s).', strjoin(defs.helmholtz_phases, ', '));
end
T = T(:);
Tc = 647.096;

switch lower(char(kind))
    case 'saturation'
        Tq = min(T, Tc - 1e-6);
        P0 = psat_aux(Tq); P0(T >= Tc) = NaN;
        fA = @(P, T_) G_rho(fluid, P, T_, 'liquid');
        fB = @(P, T_) G_rho(fluid, P, T_, 'vapor');
        [P, rA, rB] = newton_lnP(fA, fB, T, P0);
        same = ~(abs(rA ./ rB - 1) > 1e-6);
        P(same) = NaN; rA(same) = NaN; rB(same) = NaN;
    case 'sublimation'
        if ~ismember(ice, defs.solid_phases)
            error('SF_coexistence:badInput', '''%s'' is not an ice phase.', ice);
        end
        sp = sf_load_spline(fluid);
        Tmin = sp.Tc * exp(sp.knots{2}(1));
        P0 = psub_r1408(min(T, 273.16));
        fA = @(P, T_) G_rho(ice, P, T_, '');
        fB = @(P, T_) vapour(sp, fluid, P, T_, dilute, Tmin);
        [P, rA, rB] = newton_lnP(fA, fB, T, P0);
    otherwise
        error('SF_coexistence:badInput', 'kind must be ''saturation'' or ''sublimation''.');
end
out.T = T; out.P = P; out.rho_A = rA; out.rho_B = rB;
end


% -------------------------------------------------------------------------
function [G, rho] = G_rho(mat, P, T, branch)
if isempty(branch)
    o = SF_getprop([P(:) T(:)], mat, {'G', 'rho'});
else
    o = SF_getprop([P(:) T(:)], mat, {'G', 'rho'}, 'branch', branch);
end
G = o.G(:); rho = o.rho(:);
end

function [G, rho] = vapour(sp, fluid, P, T, dilute, Tmin)
persistent warned
[G, rho] = G_rho(fluid, P, T, 'vapor');
if dilute
    lo = T(:) < Tmin;
    if any(lo)
        [G(lo), rho(lo)] = sf_ideal_gas(sp, P(lo), T(lo));
        if isempty(warned)
            warned = true;             % set first: shown once even if escalated to an error
            warning('SeaFreeze:diluteExtension', ...
                ['Sublimation below %.4g K (the lowest temperature of %s) uses the dilute-vapour ' ...
                 'extension: the vapour is the surface''s ideal-gas part (Z = 1). Non-ideality ' ...
                 'there changes p_sub by ~1e-5 relative. Pass ''dilute_extension'', false for NaN ' ...
                 'instead. (Shown once per session.)'], Tmin, fluid);
        end
    end
end
end

function [P, rA, rB] = newton_lnP(fA, fB, T, P0)
lnP = log(P0(:));
live = isfinite(lnP);
rA = NaN(size(T)); rB = NaN(size(T));
for it = 1:60
    if ~any(live), break; end
    i = find(live);
    Pi = exp(lnP(i));
    [GA, rhoA] = fA(Pi, T(i));
    [GB, rhoB] = fB(Pi, T(i));
    slope = Pi * 1e6 .* (1 ./ rhoA - 1 ./ rhoB);
    step = (GA - GB) ./ slope;
    bad = ~isfinite(step);
    step(bad) = 0;
    step = max(min(step, 1), -1);
    lnP(i) = lnP(i) - step;
    lnP(i(bad)) = NaN;
    live(i) = isfinite(lnP(i)) & abs(step) > 1e-12;
end
P = exp(lnP);
ok = isfinite(P);
if any(ok)
    [~, rA(ok)] = fA(P(ok), T(ok));
    [~, rB(ok)] = fB(P(ok), T(ok));
end
end

function p = psat_aux(T)
% IAPWS-95 auxiliary saturation pressure (MPa), Wagner & Pruss 2002 eq. 2.5
Tc = 647.096; Pc = 22.064;
a = [-7.85951783 1.84408259 -11.7866497 22.6807411 -15.9618719 1.80122502];
th = 1 - T / Tc;
p = Pc * exp(Tc ./ T .* (a(1)*th + a(2)*th.^1.5 + a(3)*th.^3 + a(4)*th.^3.5 + a(5)*th.^4 + a(6)*th.^7.5));
end

function p = psub_r1408(T)
% IAPWS R14-08 sublimation pressure of ice Ih (MPa)
a = [-0.212144006e2 0.273203819e2 -0.610598130e1];
b = [0.333333333e-2 0.120666667e1 0.170333333e1];
th = T / 273.16;
s = a(1) * th.^b(1) + a(2) * th.^b(2) + a(3) * th.^b(3);
p = 611.657e-6 * exp(s ./ th);
end
