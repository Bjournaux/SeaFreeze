% run_all_tests.m  — run from Matlab/ directory
%
% Runs all script-based tests and classdef unittest suites in test/.
addpath(genpath(pwd));

nfail = 0;

% --- Script-based tests (print their own pass/fail counts) ----------------
scripts = {'test_fnGval_vs_1p0', 'test_general', 'test_input_validation', ...
           'test_selective_props', 'test_SF_PhaseLines', 'test_SF_rho2P', 'test_fnFval', 'test_phase_diagrams', ...
           'test_water_Brown2026_internals', 'test_SeaFreeze_vs_python'};
for k = 1:numel(scripts)
    fprintf('\n=== %s ===\n', scripts{k});
    try
        run(scripts{k});
    catch e
        fprintf('[ERROR] %s\n', e.message); nfail = nfail + 1;
    end
end

% --- Classdef unittest tests (auto-discovered in test/) -------------------
% matlab.unittest does not exist in GNU Octave: there the classdef suites are
% skipped and sf_verify_octave (the Octave-vs-MATLAB-reference check) runs instead.
if sf_is_octave()
    fprintf('\n=== sf_verify_octave (classdef suites need matlab.unittest: skipped under Octave) ===\n');
    try
        sf_verify_octave;
    catch e
        fprintf('[ERROR] %s\n', e.message); nfail = nfail + 1;
    end
else
    fprintf('\n=== classdef unittest suites ===\n');
    results = runtests('test');
    nfail_class = sum([results.Failed]);
    nfail = nfail + nfail_class;
    fprintf('%d passed, %d failed\n', sum([results.Passed]), nfail_class);
end

% --- Summary ---------------------------------------------------------------
fprintf('\n=============================\n');
if nfail == 0
    fprintf('ALL SUITES PASSED\n');
else
    fprintf('%d suite(s) had failures\n', nfail);
end
