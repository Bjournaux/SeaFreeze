function test_water3_internals()
% Tests for the 1.2 internals behind water_Brown2026 and the phase diagrams:
%   psi_val grid mode, the sf_load_spline cache, the SeaFreeze:longRuntime and
%   SeaFreeze:diluteExtension warnings, and SF_WPD with a Helmholtz liquid.
% Runs under MATLAB and GNU Octave.
% Baptiste Journaux - 2026
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'internal'));
addpath(fullfile(fileparts(here), 'internal', 'compat'));
addpath(fileparts(here));
np = 0; nf = 0;
sp = sf_load_spline('water_Brown2026');

% =========================================================================
% 1. psi_val 'grid' mode == scattered evaluation (all regimes)
% =========================================================================
% densities span the virial continuation (< 1e-4 kg/m3), the dilute vapour,
% the dense liquid and beyond the top knot (NaN); temperatures include one
% below the lowest knot (NaN)
rg = [1e-9; 1e-6; 1e-4; 0.5; 322; 1000; 1500; 5000; 2e4];
Tg = [200; 231; 300; 647; 1500; 3e4];
nd = struct('Fr', true, 'Frr', true);
g = psi_val(sp, rg, Tg, nd, 'grid');
[R, TT] = ndgrid(rg, Tg);
s = psi_val(sp, R(:), TT(:), nd);
e1 = max(abs(g.Fr(:) - s.Fr) ./ max(abs(s.Fr), 1e-300));
e2 = max(abs(g.Frr(:) - s.Frr) ./ max(abs(s.Frr), 1e-300));
[np,nf] = check(sprintf('psi_val grid == scattered: Fr %.1e, Frr %.1e', e1, e2), ...
    isequal(size(g.Fr), [numel(rg) numel(Tg)]) && isequal(isnan(g.Fr(:)), isnan(s.Fr)) && ...
    isequal(isnan(g.Frr(:)), isnan(s.Frr)) && e1 < 1e-10 && e2 < 1e-10, np, nf);
[np,nf] = check('psi_val grid: NaN below the lowest T knot and above the top density knot', ...
    all(isnan(g.Fr(:, 1))) && all(isnan(g.Fr(end, :))), np, nf);
try
    psi_val(sp, rg, Tg, struct('F', true), 'grid'); bad = false;
catch e
    bad = strcmp(e.identifier, 'psi_val:grid');
end
[np,nf] = check('psi_val grid refuses unsupported derivatives', bad, np, nf);

% =========================================================================
% 2. sf_load_spline cache: identical content, no second file read
% =========================================================================
clear sf_load_spline
tic; a = sf_load_spline('water_Brown2026'); t1 = toc;
tic; b = sf_load_spline('water_Brown2026'); t2 = toc;
[np,nf] = check(sprintf('sf_load_spline cache: identical struct, repeat %.1f ms vs first %.1f ms', ...
    1e3 * t2, 1e3 * t1), isequal(a, b) && t2 < t1 / 5, np, nf);
c = sf_load_spline('Ih');
[np,nf] = check('sf_load_spline cache keeps materials apart', ~isequal(c.coefs, a.coefs) && isequal(sf_load_spline('Ih'), c), np, nf);

% =========================================================================
% 3. SeaFreeze:longRuntime warning: where it should, and only there
% =========================================================================
LR = 'SeaFreeze:longRuntime';
[np,nf] = check('SF_WPD_PT warns SeaFreeze:longRuntime', ...
    issues_warning(@() SF_WPD_PT('nP', 25, 'nT', 20), LR), np, nf);
[np,nf] = check('SF_WPD_rhoT warns SeaFreeze:longRuntime', ...
    issues_warning(@() SF_WPD_rhoT('nP', 60, 'nT', 20, 'nrho', 40), LR), np, nf);
[np,nf] = check('SF_WPD with water_Bollengier2019 does not warn', ...
    ~issues_warning(@() SF_WPD('liquid', 'water_Bollengier2019', 'labels', false), LR), np, nf);

% =========================================================================
% 4. SF_WPD with the Helmholtz liquid (and its warning)
% =========================================================================
w = issues_warning(@() SF_WPD('liquid', 'water_Brown2026', 'labels', false, 'meta', 'none'), LR);
st = warning('off', LR);
f = SF_WPD('liquid', 'water_Brown2026', 'labels', false, 'meta', 'none');
warning(st);
nl = numel(findall(f, 'Type', 'line'));   % SF_WPD hides its lines from findobj
[np,nf] = check(sprintf('SF_WPD(''liquid'',''water_Brown2026'') draws the diagram (%d lines) and warns', nl), ...
    nl >= 10 && w, np, nf);
close(f);
try
    SF_WPD('liquid', 'water_Brown2018'); bad = false;
catch
    bad = true;
end
[np,nf] = check('SF_WPD rejects an unsupported liquid', bad, np, nf);

% =========================================================================
% 5. SeaFreeze:diluteExtension: once per session, only below 230 K
% =========================================================================
DE = 'SeaFreeze:diluteExtension';
clear SF_coexistence
w1 = issues_warning(@() SF_coexistence('sublimation', 200), DE);
w2 = issues_warning(@() SF_coexistence('sublimation', 190), DE);
clear SF_coexistence
w3 = issues_warning(@() SF_coexistence('sublimation', 250), DE);
w4 = issues_warning(@() SF_coexistence('sublimation', 200, 'dilute_extension', false), DE);
st = warning('off', DE);
s1 = SF_coexistence('sublimation', 200);
s4 = SF_coexistence('sublimation', 200, 'dilute_extension', false);
warning(st);
[np,nf] = check('dilute warning: first use warns, second is silent', w1 && ~w2, np, nf);
[np,nf] = check('dilute warning: silent at 250 K and when switched off', ~w3 && ~w4, np, nf);
[np,nf] = check('dilute extension on by default (finite p_sub at 200 K), NaN when off', ...
    isfinite(s1.P) && isnan(s4.P), np, nf);
clear SF_coexistence

fprintf('\n%d passed, %d failed\n', np, nf);
if nf > 0, error('test_water3_internals:failed', '%d test(s) failed', nf); end
end


function warned = issues_warning(fn, id)
% True if fn issues warning `id`.  The warning is escalated to an error for the
% call: lastwarn cannot be used, because Octave records even disabled warnings
% there (any later internal warning would overwrite the one being checked).
s0 = warning('query', id);   % restore only this id: Octave keeps a per-id
                             % 'error' state through a whole-state restore
warning('error', id);
try
    out = fn(); %#ok<NASGU>
    warned = false;
catch e
    warned = strcmp(e.identifier, id);
    if ~warned, warning(s0.state, id); rethrow(e); end
end
warning(s0.state, id);
close all
end


function [np, nf] = check(name, cond, np, nf)
if cond, fprintf('[PASS] %s\n', name); np = np + 1;
else,    fprintf('[FAIL] %s\n', name); nf = nf + 1; end
end
