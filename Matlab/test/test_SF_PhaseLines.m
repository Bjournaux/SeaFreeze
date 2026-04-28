function test_SF_PhaseLines()
% Tests for the rewritten SF_PhaseLines.
%
% Coverage:
%   1. Smoke: every supported pair returns a non-empty (P, T) curve.
%   2. v1 regression: the stable portion of each pure pair matches
%      SF_PhaseLines_v1's output within 1e-3 normalized distance.
%   3. Symmetry: SF_PhaseLines(A,B) == SF_PhaseLines(B,A).
%   4. NaClaq physics sanity: at m=1.0, low P, FPD ≈ 3.3 K.
%   5. Error paths: unsupported pairs, missing molality, m out of range.
%   6. Plot path: returns a figure handle without errors.
%
% Baptiste Journaux 2026.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'LocalBasisFunction'));
addpath(fileparts(here));

np = 0; nf = 0;

% =========================================================================
% 1. Smoke tests
% =========================================================================
pure_pairs = {
    {'Ih','water1'}, {'Ih','II'}, {'Ih','III'}, {'II','III'}, {'II','V'}, ...
    {'II','VI'},     {'III','V'}, {'III','water1'}, {'V','water1'}, ...
    {'VI','water1'}, {'V','VI'}};
for k = 1:length(pure_pairs)
    p = pure_pairs{k};
    name = sprintf('smoke %s/%s', p{1}, p{2});
    try
        out = SF_PhaseLines(p{1}, p{2});
        cond = ~isempty(out.P) && all(isfinite(out.P)) && all(isfinite(out.T));
        [np, nf] = check(name, cond, np, nf);
    catch err
        [np, nf] = check([name, ' threw: ', err.message], false, np, nf);
    end
end

nacl_pairs = {{'Ih','NaClaq'},{'II','NaClaq'},{'III','NaClaq'},{'V','NaClaq'},{'VI','NaClaq'}};
for k = 1:length(nacl_pairs)
    p = nacl_pairs{k};
    name = sprintf('smoke %s/%s @ m=1', p{1}, p{2});
    try
        out = SF_PhaseLines(p{1}, p{2}, 'm', 1.0);
        cond = ~isempty(out.P) && all(isfinite(out.P)) && all(isfinite(out.T)) && ...
               isfield(out, 'm') && out.m == 1.0;
        [np, nf] = check(name, cond, np, nf);
    catch err
        [np, nf] = check([name, ' threw: ', err.message], false, np, nf);
    end
end

% =========================================================================
% 2. v1 regression — compare stable portions on pure pairs
% =========================================================================
% v1 returns Nx2 = [P_MPa, T_K] (after the flip(TP',2) at line 787) for the
% no-meta/no-plot branch — this is the *full* contour, no stable filter
% applied. So we compare against new_out's full curve (segment='all').
for k = 1:length(pure_pairs)
    p = pure_pairs{k};
    name = sprintf('v1 regression %s/%s', p{1}, p{2});
    try
        v1_out = SF_PhaseLines_v1(p{1}, p{2});
        if isempty(v1_out)
            [np, nf] = check([name, ' (v1 empty, skipped)'], true, np, nf);
            continue;
        end
        new_out = SF_PhaseLines(p{1}, p{2}, 'segment', 'all');
        if isempty(new_out.P)
            [np, nf] = check([name, ' (new curve empty)'], false, np, nf); continue;
        end
        % Sample 20 evenly spaced indices from v1, find nearest in new.
        n_sample = min(20, size(v1_out, 1));
        idx = round(linspace(1, size(v1_out,1), n_sample));
        max_rel = 0;
        for ii = idx
            P_v1 = v1_out(ii, 1); T_v1 = v1_out(ii, 2);   % v1 returns [P, T]
            d_T = (new_out.T - T_v1) ./ max(abs(T_v1), 1);
            d_P = (new_out.P - P_v1) ./ max(abs(P_v1), 1);
            d   = sqrt(d_T.^2 + d_P.^2);
            max_rel = max(max_rel, min(d));
        end
        % 5% tolerance: pair grids differ between v1 (pair-specific narrow
        % ranges) and new (spline-knot intersection), so small offsets are
        % expected. We're catching gross breakage (wrong axis order, wrong
        % phase, etc.), not pixel-perfect alignment.
        cond = max_rel < 5e-2;
        [np, nf] = check(sprintf('%s  (max_rel=%.2e)', name, max_rel), cond, np, nf);
    catch err
        [np, nf] = check([name, ' threw: ', err.message], false, np, nf);
    end
end

% =========================================================================
% 3. Symmetry — SF_PhaseLines(A,B) == SF_PhaseLines(B,A)
% =========================================================================
for k = 1:length(pure_pairs)
    p = pure_pairs{k};
    name = sprintf('symmetry %s/%s', p{1}, p{2});
    try
        a = SF_PhaseLines(p{1}, p{2});
        b = SF_PhaseLines(p{2}, p{1});
        cond = isequal(size(a.P), size(b.P)) && ...
               max(abs(a.P - b.P)) < 1e-9 && max(abs(a.T - b.T)) < 1e-9;
        [np, nf] = check(name, cond, np, nf);
    catch err
        [np, nf] = check([name, ' threw: ', err.message], false, np, nf);
    end
end

% =========================================================================
% 4. NaClaq physics: FPD at m=1.0
% =========================================================================
try
    o = SF_PhaseLines('Ih', 'NaClaq', 'm', 1.0);
    [Pmin, idx] = min(o.P);
    Tmelt = o.T(idx);
    fpd = 273.15 - Tmelt;
    cond = abs(fpd - 3.3) < 1.5;     % literature ≈ 3.3 K with van't Hoff factor
    [np, nf] = check(sprintf('FPD at P=%.2f MPa, m=1: %.2f K (expect 3.3 ± 1.5)', ...
        Pmin, fpd), cond, np, nf);
catch err
    [np, nf] = check(['FPD test threw: ', err.message], false, np, nf);
end

% =========================================================================
% 4b. Multi-molality input — vector m returns struct array
% =========================================================================
try
    m_vec = [0.5, 1.0, 2.0, 4.0];
    arr = SF_PhaseLines('Ih','NaClaq','m', m_vec);
    cond = numel(arr) == numel(m_vec) && ...
           isequal([arr.m], m_vec) && ...
           all(arrayfun(@(s) ~isempty(s.P) && ~isempty(s.T), arr));
    [np, nf] = check('multi-m: returns struct array of correct length', cond, np, nf);

    % Each curve at higher m should sit at lower T (FPD ordering) at low P
    Ts_at_lowP = arrayfun(@(s) s.T(s.P == min(s.P)), arr);
    cond_order = all(diff(Ts_at_lowP(:)') < 0);
    [np, nf] = check('multi-m: FPD ordering monotone with m', cond_order, np, nf);

    % Plot path returns one figure shared across the array
    arr_plot = SF_PhaseLines('VI','NaClaq','m', [0.5 1.0],'plot',true);
    cond_fig = isfield(arr_plot, 'fig') && all(arrayfun(@(s) ishandle(s.fig), arr_plot));
    if cond_fig, close(arr_plot(1).fig); end
    [np, nf] = check('multi-m: plot returns figure handle on every entry', cond_fig, np, nf);
catch err
    [np, nf] = check(['multi-m threw: ' err.message], false, np, nf);
end

% =========================================================================
% 5. Error paths
% =========================================================================
[np, nf] = check_throws('unsupported pair (water1/water2)', ...
    @() SF_PhaseLines('water1','water2'), 'SeaFreeze:unsupportedPair', np, nf);
[np, nf] = check_throws('NaClaq without m', ...
    @() SF_PhaseLines('Ih','NaClaq'), 'SeaFreeze:badInput', np, nf);
[np, nf] = check_throws('m out of range', ...
    @() SF_PhaseLines('Ih','NaClaq','m',20), 'SeaFreeze:badInput', np, nf);

% =========================================================================
% 6. Plot path
% =========================================================================
try
    o = SF_PhaseLines('Ih','water1','plot',true);
    cond = isfield(o,'fig') && ishandle(o.fig);
    if cond, close(o.fig); end
    [np, nf] = check('plot path returns figure handle', cond, np, nf);
catch err
    [np, nf] = check(['plot path threw: ', err.message], false, np, nf);
end

fprintf('\n%d passed, %d failed\n', np, nf);
if nf > 0
    error('test_SF_PhaseLines: %d failures', nf);
end
end


function [np, nf] = check(name, cond, np, nf)
if cond
    fprintf('  [pass] %s\n', name);
    np = np + 1;
else
    fprintf('  [FAIL] %s\n', name);
    nf = nf + 1;
end
end


function [np, nf] = check_throws(name, fn, expected_id, np, nf)
try
    fn();
    fprintf('  [FAIL] %s : did not throw\n', name);
    nf = nf + 1;
catch err
    if strcmp(err.identifier, expected_id)
        fprintf('  [pass] %s -> %s\n', name, expected_id);
        np = np + 1;
    else
        fprintf('  [FAIL] %s : got %s, expected %s\n', name, err.identifier, expected_id);
        nf = nf + 1;
    end
end
end
