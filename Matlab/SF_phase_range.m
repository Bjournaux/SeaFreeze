function rng = SF_phase_range(material)
% SF_phase_range  Knot-domain ranges for a SeaFreeze material.
%
% Reads the .mat-bundled spline for the requested material and returns the
% knot-defined valid range in pressure, temperature, and (for compositional
% splines) molality.
%
% Usage:
%   rng = SF_phase_range('Ih')
%       rng.P = [0   400]    % MPa
%       rng.T = [1   301]    % K
%
%   rng = SF_phase_range('NaClaq')
%       rng.P = [0   5000.1] % MPa
%       rng.T = [229 501]    % K
%       rng.m = [0   7.01]   % mol/kg
%
% Materials follow SF_getprop's naming.

if ~(ischar(material) || (isstring(material) && isscalar(material)))
    error('SeaFreeze:badInput', '''material'' must be a string or character vector.');
end
material = char(material);

switch material
    case 'Ih',           file = 'SeaFreeze_Gibbs.mat';                  var = 'G_iceIh';
    case 'II',           file = 'SeaFreeze_Gibbs.mat';                  var = 'G_iceII';
    case 'III',          file = 'SeaFreeze_Gibbs.mat';                  var = 'G_iceIII';
    case 'V',            file = 'SeaFreeze_Gibbs.mat';                  var = 'G_iceV';
    case 'VI',           file = 'SeaFreeze_Gibbs.mat';                  var = 'G_iceVI';
    case 'VII_X_French', file = 'SeaFreeze_Gibbs.mat';                  var = 'G_iceVII_X_French';
    case 'water1',       file = 'SeaFreeze_Gibbs.mat';                  var = 'G_H2O_2GPa_500K';
    case 'water2',       file = 'SeaFreeze_Gibbs.mat';                  var = 'G_H2O_100GPa_10000K';
    case 'water_IAPWS95',file = 'SeaFreeze_Gibbs.mat';                  var = 'G_H2O_IAPWS';
    case 'NaClaq',       file = 'SeaFreeze_Gibbs_VII_NaCl5GPa.mat';     var = 'sp_NaCl_5GPa_500K';
    otherwise
        error('SeaFreeze:unknownMaterial', 'Unknown material ''%s''.', material);
end

S  = load(file, var);
sp = S.(var);

if ~iscell(sp.knots)
    error('SeaFreeze:badInput', 'spline.knots is not a cell — unexpected layout for material %s.', material);
end

rng.P = [sp.knots{1}(1), sp.knots{1}(end)];
rng.T = [sp.knots{2}(1), sp.knots{2}(end)];
if length(sp.knots) >= 3
    rng.m = [sp.knots{3}(1), sp.knots{3}(end)];
end
end
