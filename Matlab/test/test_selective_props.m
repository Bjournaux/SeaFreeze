function test_selective_props()
% Verify that SF_getprop returns exactly the fields the caller requested:
%   - single string  -> struct has only that field
%   - cell array     -> struct has only those fields
%   - omitted/empty  -> struct has all properties valid for the material
%   - shear/Vp/Vs    -> internal rho/Ks dependencies are stripped from output
%                       unless the caller also asked for them
%
% Baptiste Journaux - 2026

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'LocalBasisFunction'));
addpath(fileparts(here));

PT_ice    = [900 255];                % single (P,T) for ice VI
PT_water  = [200 300];                % single (P,T) for liquid water
PT_grid   = {400:10:500, 240:5:260};  % grid for ice V
PTm_nacl  = {0.1:50:200, 273:10:300, [0.1 0.5]};  % grid for NaClaq

% Property reference sets (must match SF_getprop.m):
BASE   = {'G','S','U','H','A','rho','Cp','Cv','Kt','Kp','Ks','alpha','vel','Js','gamma_Gruneisen','P','T'};
SHEAR  = {'shear','Vp','Vs'};
MIX    = {'mus','muw','f','m','xs','xw','Va','Cpa','Vm','Vw','Cpm','phi','Vex','aw'};

cases = {
    % name,                          PT,       material,   props,                      expected fields
    {'single string -> single field', PT_ice,   'VI',       'rho',                      {'rho'}}
    {'cell with 1 elem',              PT_ice,   'VI',       {'G'},                      {'G'}}
    {'multi-property cell (solid)',   PT_ice,   'VI',       {'G','rho','Cp'},           {'G','rho','Cp'}}
    {'liquid: base only',             PT_water, 'water1',   {'G','rho','Cp'},           {'G','rho','Cp'}}
    {'NaClaq base subset',            PTm_nacl, 'NaClaq',   {'rho','Cp'},               {'rho','Cp'}}
    {'NaClaq mixing only',            PTm_nacl, 'NaClaq',   {'phi','aw'},               {'phi','aw'}}
    {'NaClaq mixed base+mixing',      PTm_nacl, 'NaClaq',   {'G','rho','muw','phi'},    {'G','rho','muw','phi'}}
    {'grid input, subset',            PT_grid,  'V',        {'G','rho'},                {'G','rho'}}
    % --- shear/Vp/Vs strip-internals invariants -----------------------------
    {'Vp only -> rho/Ks stripped',    PT_ice,   'VI',       'Vp',                       {'Vp'}}
    {'Vs only -> rho/Ks stripped',    PT_ice,   'VI',       'Vs',                       {'Vs'}}
    {'shear only -> rho/Ks stripped', PT_ice,   'VI',       'shear',                    {'shear'}}
    {'{Vp,rho} -> rho kept, Ks not',  PT_ice,   'VI',       {'Vp','rho'},               {'Vp','rho'}}
    {'{Vp,Ks} -> Ks kept, rho not',   PT_ice,   'VI',       {'Vp','Ks'},                {'Vp','Ks'}}
    {'{Vp,Vs,shear,rho,Ks}',          PT_ice,   'VI',       {'Vp','Vs','shear','rho','Ks'}, {'Vp','Vs','shear','rho','Ks'}}
};

n_pass = 0; n_fail = 0;

% --- Explicit-props cases ---------------------------------------------------
for i = 1:size(cases,1)
    name = cases{i}{1};
    PT   = cases{i}{2};
    mat  = cases{i}{3};
    pr   = cases{i}{4};
    want = cases{i}{5};
    try
        out = SF_getprop(PT, mat, pr);
        got = fieldnames(out)';
        [missing, extra] = field_diff(got, want);
        if isempty(missing) && isempty(extra)
            fprintf('  [pass] %-40s  fields = {%s}\n', name, strjoin(got, ','));
            n_pass = n_pass + 1;
        else
            fprintf('  [FAIL] %-40s\n         missing: {%s}  extra: {%s}\n', ...
                    name, strjoin(missing, ','), strjoin(extra, ','));
            n_fail = n_fail + 1;
        end
    catch err
        fprintf('  [FAIL] %-40s threw: %s\n', name, err.message);
        n_fail = n_fail + 1;
    end
end

% --- "All properties" cases: no props arg / empty props ----------------------
all_cases = {
    % name,                  PT,        material,    expected superset
    {'all (water1)',         PT_water,  'water1',    BASE}
    {'all (ice VI)',         PT_ice,    'VI',        [BASE, SHEAR]}
    {'all (NaClaq)',         PTm_nacl,  'NaClaq',    [BASE, MIX]}
};
for i = 1:size(all_cases,1)
    name = all_cases{i}{1};
    PT   = all_cases{i}{2};
    mat  = all_cases{i}{3};
    want = all_cases{i}{4};
    try
        out = SF_getprop(PT, mat);  % no props arg
        got = fieldnames(out)';
        % "All" mode must return at least the expected properties (and only
        % those defined as valid for the material; no foreign fields).
        missing = setdiff(want, got);
        % Allow extras only if they're in the valid set for the material.
        if strcmp(mat,'NaClaq'),  valid = [BASE, MIX];
        elseif any(strcmp(mat,{'water1','water2','water_IAPWS95'})), valid = BASE;
        else, valid = [BASE, SHEAR]; end
        extra = setdiff(got, valid);
        if isempty(missing) && isempty(extra)
            fprintf('  [pass] %-40s  %d fields returned\n', name, numel(got));
            n_pass = n_pass + 1;
        else
            fprintf('  [FAIL] %-40s\n         missing: {%s}  unexpected: {%s}\n', ...
                    name, strjoin(missing, ','), strjoin(extra, ','));
            n_fail = n_fail + 1;
        end
    catch err
        fprintf('  [FAIL] %-40s threw: %s\n', name, err.message);
        n_fail = n_fail + 1;
    end
end

% --- Numerical sanity: subset values match all-props values ------------------
% (catches accidental coupling between selective-props code path and values)
try
    full   = SF_getprop(PT_ice, 'VI');
    subset = SF_getprop(PT_ice, 'VI', {'rho','Cp','Vp'});
    diffs = [abs(full.rho-subset.rho), abs(full.Cp-subset.Cp), abs(full.Vp-subset.Vp)];
    if all(diffs < 1e-9 * max(abs([full.rho, full.Cp, full.Vp])))
        fprintf('  [pass] subset values match all-props values\n');
        n_pass = n_pass + 1;
    else
        fprintf('  [FAIL] subset values differ from all-props: %s\n', mat2str(diffs));
        n_fail = n_fail + 1;
    end
catch err
    fprintf('  [FAIL] numeric-consistency check threw: %s\n', err.message);
    n_fail = n_fail + 1;
end

fprintf('\n%d passed, %d failed\n', n_pass, n_fail);
if n_fail > 0
    error('test_selective_props:fail', '%d cases failed', n_fail);
end
end

% ---------------------------------------------------------------------------
function [missing, extra] = field_diff(got, want)
missing = setdiff(want, got);
extra   = setdiff(got, want);
end
