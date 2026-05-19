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
% Materials follow SF_getprop's naming.

if ~(ischar(material) || (isstring(material) && isscalar(material)))
    error('SeaFreeze:badInput', '''material'' must be a string or character vector.');
end
material = char(material);

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
