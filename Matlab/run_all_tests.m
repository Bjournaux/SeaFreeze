% run_all_tests.m  — run from Matlab/ directory
%
% Runs all script-based tests and classdef unittest suites in test/.
addpath(genpath(pwd));

nfail = 0;

% --- Script-based tests (print their own pass/fail counts) ----------------
scripts = {'test_fnGval_vs_1p0', 'test_general', 'test_input_validation', ...
           'test_selective_props', 'test_SF_PhaseLines', 'test_SF_rho2P'};
for k = 1:numel(scripts)
    fprintf('\n=== %s ===\n', scripts{k});
    try
        run(scripts{k});
    catch e
        fprintf('[ERROR] %s\n', e.message); nfail = nfail + 1;
    end
end

% --- Classdef unittest tests (auto-discovered in test/) -------------------
fprintf('\n=== classdef unittest suites ===\n');
results = runtests('test');
nfail_class = sum([results.Failed]);
nfail = nfail + nfail_class;
fprintf('%d passed, %d failed\n', sum([results.Passed]), nfail_class);

% --- Summary ---------------------------------------------------------------
fprintf('\n=============================\n');
if nfail == 0
    fprintf('ALL SUITES PASSED\n');
else
    fprintf('%d suite(s) had failures\n', nfail);
end
