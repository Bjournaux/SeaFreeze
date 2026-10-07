function rng = SF_phase_range(material)
% SF_phase_range  Knot-domain ranges for a SeaFreeze material.
%
% Reads the .mat-bundled spline for the requested material and returns the
% knot-defined valid range in pressure, temperature, and (for compositional
% splines) molality.
%
% Usage:
%   rng = SF_phase_range('Ih')
%       rng.P = [0     400]    % MPa
%       rng.T = [1     301]    % K
%
%   rng = SF_phase_range('NaClaq')         % stitched LP+HP 2026 (default)
%       rng.P = [0     10001]  % MPa
%       rng.T = [229   2001]   % K
%       rng.m = [0     7.01]   % mol/kg
%
%   rng = SF_phase_range('NaClaq_LP')      % 2026 low-P spline only
%   rng = SF_phase_range('NaClaq_HP')      % 2026 high-P spline only
%   rng = SF_phase_range('NaClaq_5GPa_2024') % Brown 2024 legacy
%
%   rng = SF_phase_range('water_Brown2026')         % Helmholtz spline F(rho,T)
%       rng.rho = [rho_lo rho_hi]  % kg/m^3, knot range
%       rng.T   = [T_lo   T_hi]    % K, knot range
%       rng.P   = [P_lo   P_hi]    % MPa, from sp.Prange if stored, else the
%                                  % extent of P(rho,T) over the knot box
%                                  % where dP/drho > 0.  Not every (P,T) in
%                                  % this box is reachable; SF_getprop
%                                  % returns NaN where it is not.
%
% Materials follow SF_getprop's naming.

if ~sf_ischarlike(material)
    error('SeaFreeze:badInput', '''material'' must be a string or character vector.');
end
material = sf_material_name(char(material));   % renamed materials

% Validate material name before attempting load
defs = sf_material_defs();
known_materials = defs.known_materials;
if ~ismember(material, known_materials)
    error('SeaFreeze:unknownMaterial', 'Unknown material ''%s''.', material);
end

if strcmp(material, 'NaClaq')
    % Stitched: report the intersection domain for T/m (avoids LP extrapolation
    % artifacts at phase boundaries), full P coverage LP_lo -> HP_hi.
    spLP = sf_load_spline('NaClaq_LP');
    spHP = sf_load_spline('NaClaq_HP');
    rng.P = [spLP.knots{1}(1),  spHP.knots{1}(end)];
    rng.T = [max(spLP.knots{2}(1), spHP.knots{2}(1)), ...
             min(spLP.knots{2}(end), spHP.knots{2}(end))];
    rng.m = [max(spLP.knots{3}(1), spHP.knots{3}(1)), ...
             min(spLP.knots{3}(end), spHP.knots{3}(end))];
elseif ismember(material, defs.helmholtz_phases)
    sp = sf_load_spline(material);
    if isfield(sp, 'eos') && strcmp(sp.eos, 'psi')
        % psi surface: knots are x = ln(rho/rhoc)/3, y = ln(T/Tc); below the
        % lowest density knot the surface is continued, so rho_lo = 0.
        rng.rho = [0, sp.rhoc * exp(3 * sp.knots{1}(end))];
        rng.T   = sp.Tc * exp([sp.knots{2}(1), sp.knots{2}(end)]);
        r = sp.rhoc * exp(3 * linspace(sp.knots{1}(1), sp.knots{1}(end), 200));
    else
        rng.rho = [sp.knots{1}(1), sp.knots{1}(end)];
        if isfield(sp, 'Tc')
            rng.T = sp.Tc * exp([sp.knots{2}(1), sp.knots{2}(end)]);
        else
            rng.T = [sp.knots{2}(1), sp.knots{2}(end)];
        end
        r = linspace(rng.rho(1), rng.rho(2), 200);
    end
    if isfield(sp, 'Prange')
        rng.P = sp.Prange(:).';
    else
        T = linspace(rng.T(1), rng.T(2), 100);
        out = fnFval(sp, {r, T}, {'P', 'Kt'}, 'rhoT');
        ok  = out.Kt > 0 & isfinite(out.P);
        rng.P = [min(out.P(ok)), max(out.P(ok))];
    end
else
    sp = sf_load_spline(material);
    if ~iscell(sp.knots)
        error('SeaFreeze:badInput', ...
              'spline.knots is not a cell — unexpected layout for material %s.', material);
    end
    rng.P = [sp.knots{1}(1), sp.knots{1}(end)];
    rng.T = [sp.knots{2}(1), sp.knots{2}(end)];
    if length(sp.knots) >= 3
        rng.m = [sp.knots{3}(1), sp.knots{3}(end)];
    end
end
end
