function test_general()
% General self-contained test of the MATLAB SeaFreeze package.
% No external reference (no Python required). Exercises:
%   1. Smoke test: every material runs and returns finite values
%   2. Output shape matches input shape (grid and scatter)
%   3. Grid vs scatter consistency for the same (P,T[,m]) points
%   4. Physical sanity (density, sound speed at standard conditions)
%   5. Thermodynamic identity Cp - Cv = T*alpha^2*Kt/rho (Kt in MPa)
%   6. Reech identity Ks/Kt = Cp/Cv
%   7. NaN propagation outside parametrization bounds
%   8. SF_WhichPhase end-to-end on a small (P,T) map
%   9. SF_WhichPhase NaCl mode reproduces freezing-point depression
%
% Baptiste Journaux - 2026

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'LocalBasisFunction'));
addpath(fileparts(here));

np = 0; nf = 0;
tol_num  = 1e-10;     % numerical (grid vs scatter)
tol_ref  = 1e-4;      % 0.01% vs hardcoded Python reference values
tol_phys = 1e-2;      % loose, for thermodynamic-identity sanity (away from ref pts)

% =========================================================================
% 1. Smoke test: every material runs and returns finite values.
% =========================================================================
materials = {'Ih','II','III','V','VI','VII_X_French','water1','water2', ...
             'water_IAPWS95','NaClaq'};
sample_PT = struct( ...
    'Ih',            [50  255], ...
    'II',            [300 245], ...
    'III',           [260 245], ...
    'V',             [500 250], ...
    'VI',            [900 255], ...
    'VII_X_French',  [5000 500], ...
    'water1',        [200 300], ...
    'water2',        [200 300], ...
    'water_IAPWS95', [50  300], ...
    'NaClaq',        [200 300 1.0]);
for i = 1:length(materials)
    m = materials{i};
    try
        out = SF_getprop(sample_PT.(m), m);
        cond = isstruct(out) && isfield(out,'rho') && isfinite(out.rho);
        [np,nf] = check(sprintf('smoke: %s', m), cond, np, nf);
    catch err
        [np,nf] = check(sprintf('smoke: %s (%s)', m, err.message), false, np, nf);
    end
end

% =========================================================================
% 2. Output shape matches input shape.
% =========================================================================
P = [800 900 1000 1100]; T = [250 260 270];
out = SF_getprop({P,T}, 'VI', 'rho');
[np,nf] = check('shape: ice VI grid 4x3', isequal(size(out.rho),[4 3]), np, nf);

PTs = [800 250; 900 255; 950 260; 1000 265; 1100 270];
out = SF_getprop(PTs, 'VI', 'rho');
[np,nf] = check('shape: ice VI scatter 5x1', isequal(size(out.rho),[5 1]), np, nf);

% NaClaq grid prepends a baseline molality column (sp.cutoff) for apparent-
% property calculations, so the m dimension is length(m)+1.
P = [100 200 300]; T = [280 300]; m = [0.1 0.5 1.0 2.0];
out = SF_getprop({P,T,m}, 'NaClaq', 'rho');
[np,nf] = check('shape: NaCl grid 3x2x(4+baseline)', ...
                isequal(size(out.rho),[3 2 5]), np, nf);

% =========================================================================
% 3. Grid vs scatter consistency (pure phases only; the cross-validation
%    suite covers NaCl, where the prepended baseline column complicates a
%    direct comparison).
% =========================================================================
P = [200 500 1000]; T = [250 260 270];
g = SF_getprop({P,T}, 'VI', 'rho');
s = SF_getprop(combvec_2d(P,T), 'VI', 'rho');
rel = max(abs(g.rho(:) - s.rho(:)) ./ abs(g.rho(:)));
[np,nf] = check('grid==scatter (ice VI rho)', rel < tol_num, np, nf);

P = [0.1 100 500]; T = [280 320 360];
g = SF_getprop({P,T}, 'water1', 'rho');
s = SF_getprop(combvec_2d(P,T), 'water1', 'rho');
rel = max(abs(g.rho(:) - s.rho(:)) ./ abs(g.rho(:)));
[np,nf] = check('grid==scatter (water1 rho)', rel < tol_num, np, nf);

% =========================================================================
% 4. Physical sanity at standard conditions.
% =========================================================================
out = SF_getprop([0.1 298], 'water1');
[np,nf] = check('water1 rho ~ 997 (0.1 MPa, 298 K)',  abs(out.rho - 997)  < 5,  np, nf);
[np,nf] = check('water1 vel ~ 1497 (0.1 MPa, 298 K)', abs(out.vel - 1497) < 20, np, nf);

out = SF_getprop([0.1 268], 'Ih');
[np,nf] = check('ice Ih rho ~ 917 (0.1 MPa, 268 K)', abs(out.rho - 917)  < 5,   np, nf);
% bulk sound speed (sqrt(Ks/rho)) is ~3100; Vp (includes shear) is ~3800
[np,nf] = check('ice Ih vel ~ 3100 (0.1 MPa, 268 K)', abs(out.vel - 3100) < 200, np, nf);

% Adding salt should raise solution density above pure water.
op = SF_getprop([0.1 298],     'water1', 'rho');
on = SF_getprop([0.1 298 1.0], 'NaClaq', 'rho');
[np,nf] = check('NaCl(1m) rho > pure water rho', on.rho > op.rho + 30, np, nf);

% =========================================================================
% 5. Thermodynamic identity: Cp - Cv = T * alpha^2 * Kt / rho
% =========================================================================
out = SF_getprop([100 300], 'water1');
lhs = out.Cp - out.Cv;
rhs = 300 * out.alpha^2 * out.Kt * 1e6 / out.rho;   % Kt MPa -> Pa
[np,nf] = check('Cp-Cv identity (water1)', abs(lhs-rhs)/abs(lhs) < tol_phys, np, nf);

out = SF_getprop([900 255], 'VI');
lhs = out.Cp - out.Cv;
rhs = 255 * out.alpha^2 * out.Kt * 1e6 / out.rho;
[np,nf] = check('Cp-Cv identity (ice VI)', abs(lhs-rhs)/abs(lhs) < tol_phys, np, nf);

% =========================================================================
% 6. Reech identity: Ks / Kt = Cp / Cv
% =========================================================================
out = SF_getprop([100 300], 'water1');
[np,nf] = check('Ks/Kt = Cp/Cv (water1)', ...
    abs(out.Ks/out.Kt - out.Cp/out.Cv) / (out.Ks/out.Kt) < tol_phys, np, nf);

% =========================================================================
% 7. NaN propagation outside parametrization bounds.
% =========================================================================
% Out-of-range may be NaN or Inf depending on spline edge behavior; both are
% acceptable as "non-finite", which is what callers should treat as invalid.
out = SF_getprop([1e9 5000], 'Ih', 'rho');
[np,nf] = check('Ih extreme P -> non-finite', ~isfinite(out.rho), np, nf);

out = SF_getprop([0.1 298; 1e9 5000], 'water1', 'rho');
[np,nf] = check('water1 mixed in/out of range', ...
                isfinite(out.rho(1)) && ~isfinite(out.rho(2)), np, nf);

% =========================================================================
% 8. SF_WhichPhase end-to-end (pure water).
% =========================================================================
[np,nf] = check('WhichPhase: liquid (0.1 MPa, 300 K)', ...
                SF_WhichPhase({0.1, 300}) == 0, np, nf);
[np,nf] = check('WhichPhase: ice Ih (0.1 MPa, 250 K)', ...
                SF_WhichPhase({0.1, 250}) == 1, np, nf);
P = 0.1:50:500; T = 240:5:280;
ph = SF_WhichPhase({P, T});
[np,nf] = check(sprintf('WhichPhase grid shape %dx%d', length(P), length(T)), ...
                isequal(size(ph),[length(P) length(T)]), np, nf);

% =========================================================================
% 9. NaCl freezing-point depression.
% =========================================================================
ph_pure  = SF_WhichPhase({0.1, 268});
ph_brine = SF_WhichPhase({0.1, 268, 2.0}, 'solute', 'NaCl');
[np,nf] = check('FPD: pure ice -> brine liquid at 268 K', ...
                ph_pure == 1 && ph_brine == 0, np, nf);

ph_cold = SF_WhichPhase({0.1, 240, 2.0}, 'solute', 'NaCl');
[np,nf] = check('FPD: 2 m brine freezes at 240 K', ph_cold == 1, np, nf);

fprintf('\n%d passed, %d failed\n', np, nf);
if nf > 0
    error('test_general:fail', '%d cases failed', nf);
end
end

% ---- helpers --------------------------------------------------------------
function [np,nf] = check(name, cond, np, nf)
if cond
    fprintf('  [pass] %s\n', name);
    np = np + 1;
else
    fprintf('  [FAIL] %s\n', name);
    nf = nf + 1;
end
end

function PT = combvec_2d(P, T)
[Pg, Tg] = ndgrid(P, T);
PT = [Pg(:), Tg(:)];
end
