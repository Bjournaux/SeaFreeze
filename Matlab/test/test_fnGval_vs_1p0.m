function test_fnGval_vs_1p0()
% Regression: compare new fnGval against frozen fnGval_1p0 baseline for
% every property they share, across ices / water1 / NaClaq, in both grid
% and scatter modes.
%
% Properties unique to the new fnGval (gamma_Gruneisen, Js, F, Vw, xs, xw,
% mask) are skipped in this comparison.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here),'LocalBasisFunction'));
addpath(fileparts(here));

rtol = 1e-10;
atol = 1e-12;
np = 0; nf = 0;

% Expected-drift overrides: properties that legitimately differ from v1.0
% because of the deliberate water-molar-mass upgrade
% (nw: 1000/18.01528 -> 1000/18.015268, NIST/Tillner-Roth & Friend).
% Relative drift for any nw-linear quantity is ~6.66e-7.
rtol_overrides = struct('muw', 1e-6, 'aw', 1e-6);

% Load splines once
S_pure = load(fullfile(fileparts(here),'SeaFreeze_Gibbs.mat'));
S_nacl = load(fullfile(fileparts(here),'SeaFreeze_Gibbs_VII_NaCl5GPa.mat'));
sp_NaCl = S_nacl.sp_NaCl_5GPa_500K;
if ~isfield(sp_NaCl,'MW'), sp_NaCl.MW = 58.44e-3; end
if ~isfield(sp_NaCl,'nu'), sp_NaCl.nu = 2;        end
sp_NaCl.Go = 1;

% Pure-phase test cases: { spline, label, grid_PT, scatter_PT }
P_grid = 100:50:500; T_grid = 240:5:280;
PT_scat = [100 250; 200 260; 300 270; 400 250; 500 280];

cases_pure = {
    S_pure.G_iceIh,           'iceIh',          {P_grid, T_grid}, PT_scat
    S_pure.G_iceII,           'iceII',          {200:50:600, 200:5:260}, [200 240; 400 250; 600 230]
    S_pure.G_iceIII,          'iceIII',         {200:20:340, 230:5:265}, [220 245; 280 250; 340 260]
    S_pure.G_iceV,            'iceV',           {350:20:600, 240:5:275}, [400 250; 500 260; 600 270]
    S_pure.G_iceVI,           'iceVI',          {700:50:1500, 240:5:300}, [800 255; 1000 260; 1200 280]
    S_pure.G_H2O_2GPa_500K,   'water1',         {0.1:50:500, 250:5:350}, [100 280; 200 300; 400 320]
};

% Properties the OLD fnGval supports — these are what we compare.
common_pure_props = {'G','S','U','H','A','rho','Cp','Cv','Kt','Kp','Ks','alpha','vel'};

for k = 1:size(cases_pure,1)
    sp    = cases_pure{k,1};
    label = cases_pure{k,2};
    PT_g  = cases_pure{k,3};
    PT_s  = cases_pure{k,4};

    % Grid
    rA = fnGval_1p0(sp, PT_g, common_pure_props);
    rB = fnGval(    sp, PT_g, common_pure_props);
    [np,nf] = compare_struct(['grid ' label], rA, rB, common_pure_props, rtol, atol, rtol_overrides, np, nf);

    % Scatter
    rA = fnGval_1p0(sp, PT_s, common_pure_props);
    rB = fnGval(    sp, PT_s, common_pure_props);
    [np,nf] = compare_struct(['scat ' label], rA, rB, common_pure_props, rtol, atol, rtol_overrides, np, nf);
end

% NaCl(aq) — common props between v1.0 and new
common_nacl_props = [common_pure_props, ...
    {'mus','muw','f','m','Va','Cpa','Vm','Cpm','phi','Vex','aw'}];
P_n = 0.1:50:500; T_n = 250:10:300; m_n = [0.5 1 2 4];
PT_grid_n = {P_n, T_n, m_n};
% Build a representative scatter list for NaClaq.
[Pg,Tg,Mg] = ndgrid(P_n([1 3 5]), T_n([1 3 5]), m_n);
PT_scat_n = [Pg(:) Tg(:) Mg(:)];

rA = fnGval_1p0(sp_NaCl, PT_grid_n, common_nacl_props);
rB = fnGval(    sp_NaCl, PT_grid_n, common_nacl_props);
[np,nf] = compare_struct('grid NaClaq', rA, rB, common_nacl_props, rtol, atol, rtol_overrides, np, nf);

% Scatter NaCl: the v1.0 path is documented broken for mixing quantities,
% so we restrict the comparison to the base props on scatter input.
rA = fnGval_1p0(sp_NaCl, PT_scat_n, common_pure_props);
rB = fnGval(    sp_NaCl, PT_scat_n, common_pure_props);
[np,nf] = compare_struct('scat NaClaq (base only)', rA, rB, common_pure_props, rtol, atol, rtol_overrides, np, nf);

fprintf('\n%d passed, %d failed\n', np, nf);
if nf > 0
    error('test_fnGval_vs_1p0: %d failures', nf);
end
end

function [np,nf] = compare_struct(label, rA, rB, props, rtol, atol, overrides, np, nf)
for ii = 1:length(props)
    name = props{ii};
    if ~isfield(rA,name) || ~isfield(rB,name)
        fprintf('  [skip] %s : %s missing in one side\n', label, name);
        continue
    end
    a = rA.(name); b = rB.(name);
    if ~isequal(size(a), size(b))
        fprintf('  [FAIL] %-28s %-6s shape %s vs %s\n', label, name, mat2str(size(a)), mat2str(size(b)));
        nf = nf + 1; continue
    end
    diff = abs(a(:) - b(:));
    scale = max(abs(a(:)), abs(b(:)));
    rel = diff ./ max(scale, atol);
    if isfield(overrides, name)
        this_rtol = overrides.(name);
    else
        this_rtol = rtol;
    end
    bad = (diff > atol) & (rel > this_rtol);
    if any(bad)
        fprintf('  [FAIL] %-28s %-6s max rel %.2e  max abs %.2e  (rtol=%.0e)\n', ...
            label, name, max(rel), max(diff), this_rtol);
        nf = nf + 1;
    else
        np = np + 1;
    end
end
end
