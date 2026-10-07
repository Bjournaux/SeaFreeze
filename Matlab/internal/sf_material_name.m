function name = sf_material_name(name)
%SF_MATERIAL_NAME  Current name of a SeaFreeze material.
%
%   name = sf_material_name(name)
%
%   Renamed materials (SeaFreeze 1.2) are mapped to their new name, with a
%   SeaFreeze:deprecatedMaterial warning the first time each old name is used
%   in a session; every other name passes through unchanged.  The old names
%   keep working through SeaFreeze 1.x and are removed in 2.0.
%
%     water1            ->  water_Bollengier2019   (Bollengier et al. 2019)
%     water2            ->  water_Brown2018        (Brown 2018)
%     water3            ->  water_Brown2026        (Helmholtz psi surface, lbf-thermo 2026)
%     NaClaq_LP         ->  NaClaq_Brown2026_LP
%     NaClaq_HP         ->  NaClaq_Brown2026_HP
%     NaClaq_HP_v1..v3  ->  NaClaq_Brown2026_HP_v1..v3
%     NaClaq_5GPa_2024  ->  NaClaq_Brown2024
%
%   'NaClaq' is a permanent shortcut (no warning) for the recommended
%   NaCl(aq) model, NaClaq_Brown2026.
%
%   Suppress with: warning('off', 'SeaFreeze:deprecatedMaterial')
%   ('clear sf_material_name' shows the warnings again.)
%
% Baptiste Journaux - 2026

persistent warned
OLD = {'water1', 'water2', 'water3', 'NaClaq_LP', 'NaClaq_HP', 'NaClaq_HP_v1', ...
       'NaClaq_HP_v2', 'NaClaq_HP_v3', 'NaClaq_5GPa_2024'};
NEW = {'water_Bollengier2019', 'water_Brown2018', 'water_Brown2026', 'NaClaq_Brown2026_LP', ...
       'NaClaq_Brown2026_HP', 'NaClaq_Brown2026_HP_v1', 'NaClaq_Brown2026_HP_v2', ...
       'NaClaq_Brown2026_HP_v3', 'NaClaq_Brown2024'};

if isa(name, 'string') && isscalar(name), name = char(name); end
if ~ischar(name), return; end
if strcmp(name, 'NaClaq'), name = 'NaClaq_Brown2026'; return; end   % permanent shortcut
k = find(strcmp(OLD, name), 1);
if isempty(k), return; end

if isempty(warned), warned = {}; end
if ~any(strcmp(warned, name))
    warned{end+1} = name;            % set first: shown once even if escalated to an error
    warning('SeaFreeze:deprecatedMaterial', ...
        ['Material name ''%s'' is deprecated and will be removed in SeaFreeze 2.0; ' ...
         'use ''%s''. (Shown once per session; suppress with ' ...
         'warning(''off'',''SeaFreeze:deprecatedMaterial''))'], name, NEW{k});
end
name = NEW{k};
end
