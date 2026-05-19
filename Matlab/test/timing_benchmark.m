% timing_benchmark.m  —  measure spline load and getProp call times
% Run from Matlab/ directory via:
%   matlab -batch "cd('...SeaFreezeGenAI/Matlab'); run('test/timing_benchmark.m')"

% Add Matlab/ directory to path so SeaFreeze functions are visible
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));

fprintf('=== SeaFreeze timing benchmark (%s) ===\n\n', datestr(now));

% ---- 1. sf_load_spline cold (first call, no persistent cache) ----
materials = {'Ih', 'water1', 'VI', 'NaClaq_LP', 'NaClaq_HP', 'NaClaq_5GPa_2024'};

fprintf('--- sf_load_spline (cold, persistent cache cleared before each) ---\n');
for k = 1:numel(materials)
    mat = materials{k};
    % Clear persistent cache by clearing the function
    clear sf_load_spline
    t = tic;
    try
        sf_load_spline(mat);
        elapsed = toc(t);
        fprintf('  %-22s  %.3f s\n', mat, elapsed);
    catch ME
        fprintf('  %-22s  ERROR: %s\n', mat, ME.message);
    end
end

fprintf('\n--- sf_load_spline (warm, second call with persistent cache) ---\n');
for k = 1:numel(materials)
    mat = materials{k};
    % Prime the cache
    try
        sf_load_spline(mat);
    catch
    end
    t = tic;
    try
        sf_load_spline(mat);
        elapsed = toc(t);
        fprintf('  %-22s  %.4f s\n', mat, elapsed);
    catch ME
        fprintf('  %-22s  ERROR: %s\n', mat, ME.message);
    end
end

% ---- 2. SF_getprop single-point, cold vs warm ----
fprintf('\n--- SF_getprop single point (cold = after clear sf_load_spline) ---\n');

test_cases = {
    'VI',     300, 1000, NaN;
    'water1', 300,  100, NaN;
    'Ih',     250,  100, NaN;
    'NaClaq', 300,  500, 1.0;
};

for k = 1:size(test_cases, 1)
    mat  = test_cases{k,1};
    T    = test_cases{k,2};
    P    = test_cases{k,3};
    m    = test_cases{k,4};

    clear sf_load_spline
    t = tic;
    try
        if isnan(m)
            SF_getprop([P T], mat);
        else
            SF_getprop([P T m], mat);
        end
        elapsed = toc(t);
        fprintf('  %-10s  cold  %.3f s\n', mat, elapsed);
    catch ME
        fprintf('  %-10s  cold  ERROR: %s\n', mat, ME.message);
    end

    % warm
    t = tic;
    try
        if isnan(m)
            SF_getprop([P T], mat);
        else
            SF_getprop([P T m], mat);
        end
        elapsed = toc(t);
        fprintf('  %-10s  warm  %.4f s\n', mat, elapsed);
    catch ME
        fprintf('  %-10s  warm  ERROR: %s\n', mat, ME.message);
    end
end

fprintf('\nDone.\n');
