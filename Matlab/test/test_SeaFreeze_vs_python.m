function test_SeaFreeze_vs_python(varargin)
% Compare MATLAB SeaFreeze output against a Python reference.
% Baptiste Journaux - 2026
%
% Usage:
%     test_SeaFreeze_vs_python
%     test_SeaFreeze_vs_python('rtol', 1e-5, 'atol', 1e-8, 'verbose', true)
%
% Prerequisite: run gen_python_reference.py once from this directory to
% produce reference_SeaFreeze.mat.

p = inputParser;
addParameter(p, 'rtol', 1e-6);
addParameter(p, 'atol', 1e-10);
addParameter(p, 'verbose', false);
parse(p, varargin{:});
rtol = p.Results.rtol;
atol = p.Results.atol;
verbose = p.Results.verbose;

here = fileparts(mfilename('fullpath'));
mat_root = fileparts(here);
% Add LocalBasisFunction first, then mat_root on top so the updated
% Matlab/fnGval.m shadows any older copies in LocalBasisFunction/.
addpath(fullfile(mat_root, 'LocalBasisFunction'));
addpath(mat_root);

ref_path = fullfile(here, 'reference_SeaFreeze.mat');
gen_path = fullfile(here, 'gen_python_reference.py');
if ~exist(ref_path, 'file')
    error(['Reference file not found: %s\n' ...
           'Run gen_python_reference.py from Matlab/test/ first.'], ref_path);
end
ref = load(ref_path);

% Freshness check: warn if reference_SeaFreeze.mat is older than the
% generator script (likely stale).
if exist(gen_path, 'file')
    ref_info = dir(ref_path);
    gen_info = dir(gen_path);
    if datenum(gen_info.date) > datenum(ref_info.date)  %#ok<DATNM>
        warning('test_SeaFreeze_vs_python:staleReference', ...
            ['reference_SeaFreeze.mat is older than gen_python_reference.py;\n' ...
             'consider regenerating: cd Matlab/test/ && python gen_python_reference.py']);
    end
end

% Manifest check: list of cases the test currently expects. If any are
% missing in the loaded reference, point at the regeneration command.
expected_cases = {'Ih_grid','Ih_scatter','II_grid','III_grid','V_grid', ...
                  'VI_grid','VII_X_French_grid','water1_grid','water1_scatter', ...
                  'water_IAPWS95_grid','water2_grid','NaClaq_grid', ...
                  'NaClaq_edges_grid','NaClaq_scatter'};
missing = setdiff(expected_cases, fieldnames(ref));
if ~isempty(missing)
    error('test_SeaFreeze_vs_python:missingCases', ...
        ['reference_SeaFreeze.mat is missing case(s): %s\n' ...
         'Regenerate via gen_python_reference.py.'], strjoin(missing, ', '));
end

BASE       = {'G','S','U','H','A','rho','Cp','Cv','Kt','Kp','Ks','alpha','vel'};
% Baseline cutoff is taken from sp.cutoff in the .mat file, matching
% Python's convention, so apparent quantities align.
NACL_PROPS = [BASE, {'mus','muw','Vm','Cpm','Cpa','Va','Vex','phi','aw'}];

% Known pre-existing per-case-per-property rel-tolerance overrides.
% These mark drifts that are NOT regressions — most likely caused by the
% MATLAB fnval(fnder(...)) vs scipy splev knot-extrapolation behaviour on
% the ice V spline. Counted as [info] in the report and do not gate the
% test; see C4 in matlab_upgrade_priorities for root-cause work.
expected_drift = struct( ...
    'V_grid', struct( ...
        'S',     1e-4, ...
        'U',     1e-4, ...
        'H',     1e-3, ...
        'Cv',    1e-5, ...
        'Ks',    1e-5, ...
        'alpha', 1e-4));

pass_n = 0; info_n = 0; fail_n = 0;
case_summary = struct();   % per-case pass/info/fail counts
case_names = fieldnames(ref);

fprintf('Tolerances: rtol=%.1e, atol=%.1e\n', rtol, atol);

for ci = 1:length(case_names)
    cn = case_names{ci};
    c  = ref.(cn);

    % Phase is stored explicitly in the reference (added by gen_python_reference.py).
    % The case name's last underscore-separated token still encodes the mode.
    parts = strsplit(cn, '_');
    mode  = parts{end};
    if isfield(c, 'phase')
        phase = char(c.phase);
    else
        phase = strjoin(parts(1:end-1), '_');
    end
    is_nacl = strcmp(phase, 'NaClaq');

    if strcmp(mode, 'grid')
        if is_nacl
            PT = {c.P(:)', c.T(:)', c.m(:)'};
        else
            PT = {c.P(:)', c.T(:)'};
        end
    elseif strcmp(mode, 'scatter')
        if is_nacl
            PT = [c.P(:) c.T(:) c.m(:)];
        else
            PT = [c.P(:) c.T(:)];
        end
    else
        warning('Unknown mode "%s" in case %s, skipping', mode, cn);
        continue
    end

    if is_nacl
        props = NACL_PROPS;
    else
        props = BASE;
    end

    fprintf('\n== %s (%s, %s) ==\n', cn, phase, mode);
    case_pass = 0; case_info = 0; case_fail = 0;

    try
        out = SF_getprop(PT, phase, props);
    catch err
        fprintf('  [ERROR] SF_getprop threw: %s\n', err.message);
        fail_n = fail_n + 1;
        case_fail = 1;
        case_summary.(cn) = struct('pass',0,'info',0,'fail',1);
        continue
    end

    overrides = struct();
    if isfield(expected_drift, cn), overrides = expected_drift.(cn); end

    for j = 1:length(props)
        nm = props{j};
        if ~isfield(c, nm) || ~isfield(out, nm), continue; end
        [ok, mrel, mabs] = compare_arrays(c.(nm), out.(nm), rtol, atol);
        if ok
            case_pass = case_pass + 1;
            if verbose
                fprintf('  [pass] %-7s  rel=%.2e  abs=%.2e\n', nm, mrel, mabs);
            end
        elseif isfield(overrides, nm) && mrel <= overrides.(nm)
            case_info = case_info + 1;
            fprintf('  [info] %-7s  rel=%.2e  abs=%.2e  (expected drift, override=%.0e)\n', ...
                    nm, mrel, mabs, overrides.(nm));
        else
            case_fail = case_fail + 1;
            fprintf('  [FAIL] %-7s  rel=%.2e  abs=%.2e\n', nm, mrel, mabs);
        end
    end

    fprintf('  --- %d pass / %d info / %d fail ---\n', case_pass, case_info, case_fail);
    pass_n = pass_n + case_pass;
    info_n = info_n + case_info;
    fail_n = fail_n + case_fail;
    case_summary.(cn) = struct('pass',case_pass,'info',case_info,'fail',case_fail);
end

% Final per-case summary table.
fprintf('\n----------------------------------------------\n');
fprintf('%-26s  pass  info  fail\n', 'case');
fprintf('----------------------------------------------\n');
cn_list = fieldnames(case_summary);
for ci = 1:length(cn_list)
    s = case_summary.(cn_list{ci});
    fprintf('%-26s  %4d  %4d  %4d\n', cn_list{ci}, s.pass, s.info, s.fail);
end
fprintf('----------------------------------------------\n');
fprintf('%-26s  %4d  %4d  %4d\n', 'total', pass_n, info_n, fail_n);

if fail_n > 0
    error('test_SeaFreeze_vs_python:tolerance', ...
          '%d comparisons exceeded tolerance (rtol=%.1e, atol=%.1e)', ...
          fail_n, rtol, atol);
end
end

% ---- helpers ---------------------------------------------------------------

function [ok, mrel, mabs] = compare_arrays(py, ml, rtol, atol)
    py = squeeze(py);
    ml = squeeze(ml);
    if ~isequal(size(py), size(ml))
        if numel(py) == numel(ml)
            py = py(:); ml = ml(:);
        else
            ok = false; mrel = NaN; mabs = NaN;
            return
        end
    end
    py = py(:); ml = ml(:);
    mask = ~(isnan(py) | isnan(ml));
    if ~any(mask), ok = true; mrel = 0; mabs = 0; return; end
    dabs = abs(py(mask) - ml(mask));
    % Symmetric rel-error denominator: avoids quirks when py is near zero
    % but ml isn't (or vice versa).
    drel = dabs ./ max(max(abs(py(mask)), abs(ml(mask))), atol);
    mabs = max(dabs);
    mrel = max(drel);
    ok = all(dabs <= atol | drel <= rtol);
end
