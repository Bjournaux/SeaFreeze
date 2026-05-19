function sp = sf_load_spline(material)
% SF_LOAD_SPLINE  Load a SeaFreeze Gibbs spline by material code.
%
% Returns the spline struct for the requested material from the per-spline
% folder under Matlab/splines/.  All files use the standardised variable name
% `sp` so callers do not need to know the original source variable names.
%
% Usage:
%   sp = sf_load_spline('Ih')             % ice Ih
%   sp = sf_load_spline('II')             % ice II
%   sp = sf_load_spline('III')            % ice III
%   sp = sf_load_spline('V')              % ice V
%   sp = sf_load_spline('VI')             % ice VI
%   sp = sf_load_spline('VII_X_French')   % ice VII/X (French & Redmer 2015)
%   sp = sf_load_spline('water1')         % Bollengier 2019 (<=2 GPa, <=500 K)
%   sp = sf_load_spline('water2')         % Brown 2018 (up to 100 GPa)
%   sp = sf_load_spline('water_IAPWS95')  % IAPWS95, Wagner & Pruss 2002
%   sp = sf_load_spline('NaClaq_LP')      % 2026 low-P  NaCl(aq) LBF spline
%   sp = sf_load_spline('NaClaq_HP')      % 2026 high-P NaCl(aq) LBF spline (r3)
%   sp = sf_load_spline('NaClaq_5GPa_2024') % Brown 2024 NaCl(aq) spline (legacy)
%   sp = sf_load_spline('NaClaq_HP_v1')  % March2026 HP alt. fit v1
%   sp = sf_load_spline('NaClaq_HP_v2')  % March2026 HP alt. fit v2
%   sp = sf_load_spline('NaClaq_HP_v3')  % March2026 HP alt. fit v3
%
% Note: the stitched LP+HP default ('NaClaq' in SF_getprop) is not a single
% file and cannot be loaded here — use SF_getprop or SF_NaCl_stitch directly.
%
% The loader uses a persistent lookup table so repeated calls for the same
% material hit the MATLAB file cache rather than re-parsing the map.
%
% See also: SF_getprop, SF_phase_range, SF_WhichPhase
%
% Baptiste Journaux — 2026

persistent MAP
if isempty(MAP)
    MAP = build_map();
end

if ~(ischar(material) || (isstring(material) && isscalar(material)))
    error('sf_load_spline:badInput', ...
          '''material'' must be a character vector or scalar string.');
end
material = char(material);

if ~isKey(MAP, material)
    known = strjoin(keys(MAP), ', ');
    error('sf_load_spline:unknown', ...
          'Unknown material ''%s''.\nKnown materials: %s', material, known);
end

entry = MAP(material);           % {subfolder, filename}
spdir = fullfile(fileparts(mfilename('fullpath')), '..', 'splines');
S  = load(fullfile(spdir, entry{1}, entry{2}), 'sp');
sp = S.sp;

end  % main function


% -------------------------------------------------------------------------
function MAP = build_map()
% Returns a containers.Map: material_code → {subfolder, matfile}

keys_list = {'Ih', 'II', 'III', 'V', 'VI', 'VII_X_French', ...
             'water1', 'water2', 'water_IAPWS95', ...
             'NaClaq_LP', 'NaClaq_HP', 'NaClaq_5GPa_2024', ...
             'NaClaq_HP_v1', 'NaClaq_HP_v2', 'NaClaq_HP_v3'};

vals_list = {
    {'ice_Ih',              'ice_Ih.mat'}
    {'ice_II',              'ice_II.mat'}
    {'ice_III',             'ice_III.mat'}
    {'ice_V',               'ice_V.mat'}
    {'ice_VI',              'ice_VI.mat'}
    {'ice_VII_X_French',    'ice_VII_X_French.mat'}
    {'water_Bollengier',    'water_Bollengier.mat'}
    {'water_Brown',         'water_Brown.mat'}
    {'water_IAPWS95',       'water_IAPWS95.mat'}
    {'NaCl_aq_LP_2026',     'NaCl_aq_LP_2026.mat'}
    {'NaCl_aq_HP_2026',     'NaCl_aq_HP_2026.mat'}
    {'NaCl_aq_Brown2024',   'NaCl_aq_Brown2024.mat'}
    {'NaCl_aq_HP_2026_v1',  'NaCl_aq_HP_2026_v1.mat'}
    {'NaCl_aq_HP_2026_v2',  'NaCl_aq_HP_2026_v2.mat'}
    {'NaCl_aq_HP_2026_v3',  'NaCl_aq_HP_2026_v3.mat'}
};

MAP = containers.Map(keys_list, vals_list);
end
