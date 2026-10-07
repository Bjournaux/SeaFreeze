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
%     water1  ->  water_Bollengier2019   (Bollengier et al. 2019)
%     water2  ->  water_Brown2018        (Brown 2018)
%     water3  ->  water_Brown2026        (Helmholtz psi surface, lbf-thermo 2026)
%
%   Suppress with: warning('off', 'SeaFreeze:deprecatedMaterial')
%   ('clear sf_material_name' shows the warnings again.)
%
% Baptiste Journaux - 2026

persistent warned
OLD = {'water1', 'water2', 'water3'};
NEW = {'water_Bollengier2019', 'water_Brown2018', 'water_Brown2026'};

if isa(name, 'string') && isscalar(name), name = char(name); end
if ~ischar(name), return; end
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
