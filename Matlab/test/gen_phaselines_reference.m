function gen_phaselines_reference(varargin)
% Generate a .mat reference file from the MATLAB SF_PhaseLines implementation,
% for consumption by Python/seafreeze/test/test_phaselines_vs_matlab.py.
%
% Baptiste Journaux - 2026
%
% Run once from this directory (Matlab/test/):
%
%     gen_phaselines_reference
%
% Saves reference_phaselines.mat with one struct per pair, each containing
% the (P, T) of the equilibrium curve, the stable mask, and (for NaClaq
% pairs) the molality. Both implementations evaluate on the *same* dense
% T grid (0.005 K) so that linear-interpolation error in the contour
% extraction is far below the 0.01 K cross-validation tolerance.
%
% Optional arguments:
%     gen_phaselines_reference('dT', 0.005)    % T-grid spacing (K)

p = inputParser;
addParameter(p, 'dT', 0.005, @(x) isnumeric(x) && isscalar(x) && x > 0);
parse(p, varargin{:});
dT = p.Results.dT;

here = fileparts(mfilename('fullpath'));
mat_root = fileparts(here);
addpath(fullfile(mat_root, 'LocalBasisFunction'));
addpath(mat_root);

% --- Pure-water pairs (all 12, no VII_X_French) -----------------------------
PURE_PAIRS = {
    {'Ih',  'water1'};
    {'Ih',  'II'};
    {'Ih',  'III'};
    {'II',  'III'};
    {'II',  'V'};
    {'II',  'VI'};
    {'II',  'water1'};
    {'III', 'V'};
    {'III', 'water1'};
    {'V',   'water1'};
    {'V',   'VI'};
    {'VI',  'water1'};
};

% --- NaClaq cases (3 ices x 2 molalities) -----------------------------------
NACL_PAIRS = {
    {'Ih',  'NaClaq', 0.5};
    {'Ih',  'NaClaq', 2.0};
    {'III', 'NaClaq', 0.5};
    {'III', 'NaClaq', 2.0};
    {'V',   'NaClaq', 0.5};
    {'V',   'NaClaq', 2.0};
};

cases = struct();

fprintf('Generating SF_PhaseLines reference (dT=%.4f K)...\n', dT);

for k = 1:length(PURE_PAIRS)
    matA = PURE_PAIRS{k}{1};
    matB = PURE_PAIRS{k}{2};
    label = sprintf('%s_%s', matA, matB);

    % Use the spline-domain intersection T range with the requested step.
    rA = SF_phase_range(matA);
    rB = SF_phase_range(matB);
    T_lo = max([rA.T(1), rB.T(1), 1.0]);
    T_hi = min(rA.T(2), rB.T(2));
    T_grid = T_lo:dT:T_hi;

    fprintf('  %-22s T=[%.2f,%.2f]  nT=%d ...', label, T_lo, T_hi, numel(T_grid));
    t0 = tic;
    out = SF_PhaseLines(matA, matB, 'T', T_grid, 'segment', 'all');
    fprintf(' %d pts in %.2fs\n', numel(out.P), toc(t0));

    cases.(label) = struct( ...
        'matA',   matA, ...
        'matB',   matB, ...
        'P',      out.P(:), ...
        'T',      out.T(:), ...
        'stable', out.stable(:), ...
        'm',      NaN, ...
        'dT',     dT);
end

for k = 1:length(NACL_PAIRS)
    matA = NACL_PAIRS{k}{1};
    matB = NACL_PAIRS{k}{2};
    m_val = NACL_PAIRS{k}{3};
    label = sprintf('%s_%s_m%s', matA, matB, ...
        strrep(sprintf('%g', m_val), '.', 'p'));

    rA = SF_phase_range(matA);
    rB = SF_phase_range(matB);
    T_lo = max([rA.T(1), rB.T(1), 1.0]);
    T_hi = min(rA.T(2), rB.T(2));
    T_grid = T_lo:dT:T_hi;

    fprintf('  %-22s m=%-4g  T=[%.2f,%.2f]  nT=%d ...', ...
        label, m_val, T_lo, T_hi, numel(T_grid));
    t0 = tic;
    out = SF_PhaseLines(matA, matB, 'm', m_val, 'T', T_grid, 'segment', 'all');
    fprintf(' %d pts in %.2fs\n', numel(out.P), toc(t0));

    cases.(label) = struct( ...
        'matA',   matA, ...
        'matB',   matB, ...
        'P',      out.P(:), ...
        'T',      out.T(:), ...
        'stable', out.stable(:), ...
        'm',      m_val, ...
        'dT',     dT);
end

out_path = fullfile(here, 'reference_phaselines.mat');
save(out_path, '-struct', 'cases');
n = numel(fieldnames(cases));
fprintf('\nWrote %d cases to %s\n', n, out_path);
end
