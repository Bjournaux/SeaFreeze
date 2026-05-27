function test_SF_rho2P()
% Tests for SF_rho2P — pressure-from-density inversion.
% Baptiste Journaux - 2026
%
% Strategy: round-trip consistency.  For each material, draw a set of
% (P,T[,m]) reference points inside the spline domain, compute rho via
% SF_getprop, invert back to P with SF_rho2P, and verify the residual is
% within tolerance (0.05 MPa for most cases).
%
% Additional tests: output shape, NaN propagation, bad-input errors.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'internal'));
addpath(fileparts(here));

np = 0; nf = 0;
tol_P = 0.05;    % MPa — acceptable round-trip error

% =========================================================================
% 1. Round-trip: water1 scatter
% =========================================================================
P_ref = [0.1; 100; 500; 1000; 2000];
T_ref = [280; 300; 320; 340; 355];
out = SF_getprop([P_ref T_ref], 'water1', {'rho','Kt'});
P_rec = SF_rho2P(out.rho, T_ref, 'water1');
[np,nf] = check('roundtrip: water1 scatter (5 pts)', ...
    isequal(size(P_rec), size(P_ref)) && all(isfinite(P_rec)) && max(abs(P_rec - P_ref)) < tol_P, np, nf);

% =========================================================================
% 2. Round-trip: ice VI scatter
% =========================================================================
P_ref = [700; 900; 1100; 1500; 2000];
T_ref = [250; 255; 260; 270; 280];
out = SF_getprop([P_ref T_ref], 'VI', 'rho');
P_rec = SF_rho2P(out.rho, T_ref, 'VI');
[np,nf] = check('roundtrip: ice VI scatter (5 pts)', ...
    isequal(size(P_rec), size(P_ref)) && all(isfinite(P_rec)) && max(abs(P_rec - P_ref)) < tol_P, np, nf);

% =========================================================================
% 3. Round-trip: ice Ih scatter
% =========================================================================
P_ref = [0.1; 50; 100; 200; 350];
T_ref = [260; 265; 268; 270; 271];
out = SF_getprop([P_ref T_ref], 'Ih', 'rho');
P_rec = SF_rho2P(out.rho, T_ref, 'Ih');
[np,nf] = check('roundtrip: ice Ih scatter (5 pts)', ...
    isequal(size(P_rec), size(P_ref)) && all(isfinite(P_rec)) && max(abs(P_rec - P_ref)) < tol_P, np, nf);

% =========================================================================
% 4. Round-trip: NaClaq scatter at m = 1 mol/kg
% =========================================================================
P_ref = [0.1; 200; 500; 1500; 5000];
T_ref = [280; 300; 320; 350;  400];
m_ref = 1.0;
out = SF_getprop([P_ref T_ref repmat(m_ref, 5, 1)], 'NaClaq', 'rho');
P_rec = SF_rho2P(out.rho, T_ref, 'NaClaq', m_ref);
[np,nf] = check('roundtrip: NaClaq scatter m=1 mol/kg', ...
    isequal(size(P_rec), size(P_ref)) && all(isfinite(P_rec)) && max(abs(P_rec - P_ref)) < tol_P, np, nf);

% =========================================================================
% 5. Round-trip: NaClaq_5GPa_2024 legacy spline
% =========================================================================
P_ref = [100; 500; 1000; 2000; 3000];
T_ref = [280; 300; 320; 350; 400];
m_ref = 2.0;
out = SF_getprop([P_ref T_ref repmat(m_ref, 5, 1)], 'NaClaq_5GPa_2024', 'rho');
P_rec = SF_rho2P(out.rho, T_ref, 'NaClaq_5GPa_2024', m_ref);
[np,nf] = check('roundtrip: NaClaq_5GPa_2024 scatter m=2 mol/kg', ...
    isequal(size(P_rec), size(P_ref)) && all(isfinite(P_rec)) && max(abs(P_rec - P_ref)) < tol_P, np, nf);

% =========================================================================
% 6. Output shape: row vector input -> same row vector output
% =========================================================================
P_row = [100 500 1000 1500 2000];
T_row = [280 300 320  340  355];
out = SF_getprop([P_row(:) T_row(:)], 'water1', 'rho');
rho_row = reshape(out.rho, 1, []);   % force row
P_rec = SF_rho2P(rho_row, T_row, 'water1');
[np,nf] = check('output shape: row vector preserved', ...
    isequal(size(P_rec), [1 5]), np, nf);

% =========================================================================
% 7. Output shape: scalar inputs -> scalar output
% =========================================================================
out = SF_getprop([200 300], 'water1', 'rho');
P_sc = SF_rho2P(out.rho, 300, 'water1');
[np,nf] = check('output shape: scalar in -> scalar out', isscalar(P_sc), np, nf);

% =========================================================================
% 8. Scalar T broadcast against vector rho
% =========================================================================
P_ref = [100; 300; 600; 900; 1200];
T_fixed = 300;
out = SF_getprop([P_ref repmat(T_fixed, 5, 1)], 'water1', 'rho');
P_rec = SF_rho2P(out.rho, T_fixed, 'water1');
[np,nf] = check('scalar T broadcast: water1 (5 pts)', ...
    all(isfinite(P_rec)) && max(abs(P_rec - P_ref)) < tol_P, np, nf);

% =========================================================================
% 9. NaN for rho_target below domain minimum (less dense than allowed)
% =========================================================================
rho_lo = SF_getprop([0.1 373], 'water1', 'rho');   % hot water at low P ~ 958 kg/m³
P_test = SF_rho2P(rho_lo.rho - 200, 373, 'water1');  % target much less dense
[np,nf] = check('NaN for rho below domain min', isnan(P_test), np, nf);

% =========================================================================
% 10. NaN for rho_target above domain maximum (denser than allowed)
% =========================================================================
rho_hi = SF_getprop([2300 200], 'water1', 'rho');   % high-P cold water
P_test = SF_rho2P(rho_hi.rho + 200, 200, 'water1');
[np,nf] = check('NaN for rho above domain max', isnan(P_test), np, nf);

% =========================================================================
% 11. Physics sanity: higher rho -> higher P (isothermal compression)
% =========================================================================
T_fix = 300;
P_grid = [100; 500; 1000; 1500];
out = SF_getprop([P_grid repmat(T_fix, 4, 1)], 'water1', 'rho');
% rho increases with P; SF_rho2P must recover in the correct order
P_rec = SF_rho2P(out.rho, T_fix, 'water1');
[np,nf] = check('physics: higher rho -> higher P (water1)', ...
    all(diff(P_rec) > 0), np, nf);

% =========================================================================
% 12. Custom initial guess (P0) keyword
% =========================================================================
P_ref = [800; 900; 1000];
T_ref = [255; 258; 262];
out = SF_getprop([P_ref T_ref], 'VI', 'rho');
P_rec = SF_rho2P(out.rho, T_ref, 'VI', 'P0', 1000);
[np,nf] = check('custom P0 keyword: ice VI (3 pts)', ...
    all(isfinite(P_rec)) && max(abs(P_rec - P_ref)) < tol_P, np, nf);

% =========================================================================
% 13. Custom tolerance keyword
% =========================================================================
out = SF_getprop([500 300], 'water1', 'rho');
P_tight = SF_rho2P(out.rho, 300, 'water1', 'tol', 1e-4);
[np,nf] = check('tight tol (1e-4 MPa): water1 single pt', ...
    abs(P_tight - 500) < 0.001, np, nf);

% =========================================================================
% 14. Bad input: unknown material
% =========================================================================
try
    SF_rho2P(1000, 300, 'IceX');
    [np,nf] = check('error: unknown material', false, np, nf);
catch err
    [np,nf] = check('error: unknown material', ...
        strcmp(err.identifier,'SF_rho2P:unknownMaterial'), np, nf);
end

% =========================================================================
% 15. Bad input: NaClaq without m
% =========================================================================
try
    SF_rho2P(1000, 300, 'NaClaq');
    [np,nf] = check('error: NaClaq no m', false, np, nf);
catch err
    [np,nf] = check('error: NaClaq no m', ...
        strcmp(err.identifier,'SF_rho2P:badInput'), np, nf);
end

% =========================================================================
% 16. Round-trip: water2 (high-T region where spline is monotone)
% =========================================================================
% water2 (Brown 2018) targets extreme conditions; it is non-monotone in P at
% low-moderate T (<500 K), so the inversion is only reliable at high T.
% Use T=1000 K, P=[200,1000] MPa where density increases monotonically.
P_ref = [200; 400; 600; 800; 1000];
T_ref = repmat(1000, 5, 1);
out = SF_getprop([P_ref T_ref], 'water2', 'rho');
P_rec = SF_rho2P(out.rho, T_ref, 'water2');
[np,nf] = check('roundtrip: water2 (T=1000 K, 5 pts)', ...
    isequal(size(P_rec), size(P_ref)) && all(isfinite(P_rec)) && max(abs(P_rec - P_ref)) < 1.0, np, nf);

% =========================================================================
% 17. Low-P regression: ice Ih at 0.1 MPa (1 bar) — original bug report
% =========================================================================
% Before the fix, P_search_lo = max(0.1, 400*0.001) = 0.4 MPa for ice Ih,
% so any target P < 0.4 MPa was silently returned as NaN.
P_ref = [0.1; 0.5; 1.0; 5.0; 10.0];
T_ref = repmat(268, 5, 1);
out = SF_getprop([P_ref T_ref], 'Ih', 'rho');
P_rec = SF_rho2P(out.rho, T_ref, 'Ih');
[np,nf] = check('roundtrip: ice Ih low-P (0.1-10 MPa)', ...
    all(isfinite(P_rec)) && max(abs(P_rec - P_ref)) < tol_P, np, nf);

% =========================================================================
fprintf('\n%d passed, %d failed\n', np, nf);
if nf > 0
    error('test_SF_rho2P:fail', '%d test(s) failed', nf);
end
end

% ---- helper ---------------------------------------------------------------
function [np,nf] = check(name, cond, np, nf)
if cond
    fprintf('  [pass] %s\n', name);
    np = np + 1;
else
    fprintf('  [FAIL] %s\n', name);
    nf = nf + 1;
end
end
