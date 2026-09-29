% gen_water3_reference.m — reference data for the Helmholtz water3 tests.
%
% Writes, in test/fixtures/:
%   psi_reference.mat           psiH2O_val.m (lbf-thermo) at scattered, gridded
%                               and edge (rho,T) states — the independent
%                               reference for psi_val / fnFval and the Python
%                               port.  Needs lbf-thermo's psiH2O_val.m and the
%                               Curve Fitting Toolbox (fnval/fnder); skipped
%                               otherwise.  Set LBF_PSI_DIR below.
%   water3_getprop_reference.mat SF_getprop(...,'water3') in (P,T) grid and
%                               scatter mode — the MATLAB side of the Python
%                               parity test (test_helmholtz.py).
%
% Run from Matlab/:  run('test/gen_water3_reference.m')

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root); addpath(fullfile(root, 'internal')); addpath(fullfile(root, 'internal', 'compat'));
outdir = fullfile(here, 'fixtures');
if ~exist(outdir, 'dir'), mkdir(outdir); end

% ---- 1. psiH2O_val reference ------------------------------------------------
LBF_PSI_DIR = fullfile(getenv('HOME'), 'Documents', 'GitHub', 'lbf-thermo', ...
                       'materials', 'H2O_liquid', 'fitting', 'psi_2026');
if exist(fullfile(LBF_PSI_DIR, 'psiH2O_val.m'), 'file') && ~isempty(which('fnval'))
    addpath(LBF_PSI_DIR);
    sp = sf_load_spline('water3');
    rng(1);
    n = 400;
    rho = [900 + 500*rand(n,1); 1e-3 + 10*rand(60,1); 300 + 500*rand(60,1); ...
           1e-5*10.^(rand(20,1)); 2000 + 4000*rand(40,1)];
    T   = [235 + 400*rand(n,1); 300 + 1500*rand(60,1); 600 + 800*rand(60,1); ...
           250 + 400*rand(20,1); 1000 + 8000*rand(40,1)];
    ref.rho = rho; ref.T = T; ref.out = psiH2O_val(sp, [rho T]);
    rg = linspace(950, 1250, 31); Tg = linspace(240, 500, 27);
    ref.grid.rho = rg; ref.grid.T = Tg; ref.grid.out = psiH2O_val(sp, {rg, Tg});
    ref.edge.rho = [1e-6; 20000; 1000; 1000]; ref.edge.T = [300; 300; 200; 200000];
    ref.edge.out = psiH2O_val(sp, [ref.edge.rho ref.edge.T]);
    save(fullfile(outdir, 'psi_reference.mat'), 'ref', '-v7');
    fprintf('Wrote psi_reference.mat\n');
else
    fprintf('psiH2O_val.m / Curve Fitting Toolbox not available: psi_reference.mat not regenerated\n');
end

% ---- 2. SF_getprop water3 (P,T) reference -----------------------------------
P = [0.1 1 10 50 100 200 500 1000 2000 5000];
T = [240 260 273.16 300 350 400 500 700 1000];
g = SF_getprop({P, T}, 'water3');
[Pm, Tm] = ndgrid(P, T);
s = SF_getprop([Pm(:) Tm(:)], 'water3');
w3.P = P; w3.T = T; w3.grid = g; w3.scatter = s;
w3.sat = SF_coexistence('saturation', [273.16 300 400 500 600 640]);
w3.sub = SF_coexistence('sublimation', [180 200 230 250 273.16], 'dilute_extension', true);
rT = SF_getprop({[1000 1050 1100 1200], [280 300 330 360]}, 'water1', {'P','G','Cp'}, 'input', 'rhoT');
w3.rhoT_water1 = struct('rho', [1000 1050 1100 1200], 'T', [280 300 330 360], 'P', rT.P, 'G', rT.G, 'Cp', rT.Cp);
save(fullfile(outdir, 'water3_getprop_reference.mat'), 'w3', '-v7');
fprintf('Wrote water3_getprop_reference.mat\n');
