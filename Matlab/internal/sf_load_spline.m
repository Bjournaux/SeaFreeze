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
% Octave: ten of the spline files are stored as MATLAB v7.3 (HDF5), in which
% sp.knots is a cell array held behind HDF5 object references that Octave's
% `load` cannot follow.  Under Octave this function therefore looks in
% Matlab/splines_octave/ first, where those ten files are mirrored in MAT v7
% format (see tools/convert_splines_for_octave.py), and falls back to
% Matlab/splines/ for the rest, which are already v7 and load natively.
%
% The lookup table is built once into a persistent variable so repeated calls
% for the same material hit the MATLAB/Octave file cache rather than
% re-building the map.
%
% See also: SF_getprop, SF_phase_range, SF_WhichPhase, sf_is_octave
%
% Baptiste Journaux — 2026

persistent MAP
if isempty(MAP)
    MAP = build_map();
end

if ~sf_ischarlike(material)
    error('sf_load_spline:badInput', ...
          '''material'' must be a character vector or scalar string.');
end
material = char(material);

row = find(strcmp(MAP(:,1), material), 1);
if isempty(row)
    known = strjoin(MAP(:,1).', ', ');
    error('sf_load_spline:unknown', ...
          'Unknown material ''%s''.\nKnown materials: %s', material, known);
end

subfolder = MAP{row, 2};
filename  = MAP{row, 3};

root = fileparts(fileparts(mfilename('fullpath')));   % .../Matlab

% Under Octave prefer the converted copy when one exists for this spline.
matfile = '';
if sf_is_octave()
    candidate = fullfile(root, 'splines_octave', subfolder, filename);
    if exist(candidate, 'file') == 2
        matfile = candidate;
    end
end
if isempty(matfile)
    matfile = fullfile(root, 'splines', subfolder, filename);
end

S  = load(matfile, 'sp');
sp = S.sp;

end  % main function


% -------------------------------------------------------------------------
function MAP = build_map()
% Returns an N-by-3 cell table: {material_code, subfolder, matfile}.
%
% A plain cell table rather than containers.Map: the lookup is a strcmp over
% fifteen entries either way, and Octave's containers.Map is a classdef
% reimplementation with its own quirks that this does not need to depend on.

MAP = {
    'Ih',               'ice_Ih',              'ice_Ih.mat'
    'II',               'ice_II',              'ice_II.mat'
    'III',              'ice_III',             'ice_III.mat'
    'V',                'ice_V',               'ice_V.mat'
    'VI',               'ice_VI',              'ice_VI.mat'
    'VII_X_French',     'ice_VII_X_French',    'ice_VII_X_French.mat'
    'water1',           'water_Bollengier',    'water_Bollengier.mat'
    'water2',           'water_Brown',         'water_Brown.mat'
    'water_IAPWS95',    'water_IAPWS95',       'water_IAPWS95.mat'
    'NaClaq_LP',        'NaCl_aq_LP_2026',     'NaCl_aq_LP_2026.mat'
    'NaClaq_HP',        'NaCl_aq_HP_2026',     'NaCl_aq_HP_2026.mat'
    'NaClaq_5GPa_2024', 'NaCl_aq_Brown2024',   'NaCl_aq_Brown2024.mat'
    'NaClaq_HP_v1',     'NaCl_aq_HP_2026_v1',  'NaCl_aq_HP_2026_v1.mat'
    'NaClaq_HP_v2',     'NaCl_aq_HP_2026_v2',  'NaCl_aq_HP_2026_v2.mat'
    'NaClaq_HP_v3',     'NaCl_aq_HP_2026_v3',  'NaCl_aq_HP_2026_v3.mat'
};
end
