function gen_getProp_reference()
% Generate a .mat reference file from the MATLAB SF_getprop implementation,
% for consumption by Python/seafreeze/test/test_getProp_vs_matlab.py.
%
% Baptiste Journaux - 2026
%
% Run once from this directory (Matlab/test/):
%
%     gen_getProp_reference
%
% Saves reference_getProp.mat with one struct per case.  Each struct
% contains all SF_getprop output fields plus the input PT(m) arrays so
% the Python test can reproduce the exact same call without hard-coding.
%
% Cases
% -----
%   VI_scatter       ice VI,   scatter  5 pts
%   VI_grid          ice VI,   grid     5×4
%   Ih_scatter       ice Ih,   scatter  3 pts
%   water1_scatter   water1,   scatter  5 pts
%   water1_grid      water1,   grid     5×4
%   NaClaq_scatter   NaClaq_5GPa_2024,  scatter  4 pts  (P,T,m)
%   NaClaq_grid      NaClaq_5GPa_2024,  grid     4×4×3  {P,T,m}
%   water3_scatter   water3 (Helmholtz), scatter 5 pts (vapour, liquid, supercritical)
%   water3_grid      water3,   grid     4×4
%   water3_rhoT      water3,   (rho,T) scatter 4 pts ('input','rhoT')
%   water1_rhoT      water1,   (rho,T) grid 3×3 (Gibbs spline inverted by SF_rho2P)

here = fileparts(mfilename('fullpath'));
mat_root = fileparts(here);
addpath(fullfile(mat_root, 'internal'));
addpath(fullfile(mat_root, 'internal', 'compat'));
addpath(mat_root);

cases = struct();

fprintf('Generating SF_getprop reference data...\n');

% ---- ice VI scatter ----------------------------------------------------------
P_vec = [800 900 1000 1100 1200]';
T_vec = [250 255  260  265  270]';
pts   = [P_vec, T_vec];   % N×2 scatter matrix
out   = SF_getprop(pts, 'VI');
cases.VI_scatter = add_pt_input(out, P_vec, T_vec, [], 'VI');
fprintf('  VI_scatter          (%d pts)\n', numel(P_vec));

% ---- ice VI grid -------------------------------------------------------------
P_grid = 800:100:1200;          % 1×5
T_grid = 250:5:265;             % 1×4
out    = SF_getprop({P_grid, T_grid}, 'VI');
cases.VI_grid = add_grid_input(out, P_grid, T_grid, [], 'VI');
fprintf('  VI_grid             (%d×%d)\n', numel(P_grid), numel(T_grid));

% ---- ice Ih scatter ----------------------------------------------------------
P_vec = [1 10 50]';
T_vec = [250 260 270]';
pts   = [P_vec, T_vec];
out   = SF_getprop(pts, 'Ih');
cases.Ih_scatter = add_pt_input(out, P_vec, T_vec, [], 'Ih');
fprintf('  Ih_scatter          (%d pts)\n', numel(P_vec));

% ---- water1 scatter ----------------------------------------------------------
P_vec = [100 200 500 1000 2000]';
T_vec = [280 300 320 340  350]';
pts   = [P_vec, T_vec];
out   = SF_getprop(pts, 'water1');
cases.water1_scatter = add_pt_input(out, P_vec, T_vec, [], 'water1');
fprintf('  water1_scatter      (%d pts)\n', numel(P_vec));

% ---- water1 grid -------------------------------------------------------------
P_grid = [100 500 1000 2000 2200];
T_grid = [280 300 320  340];
out    = SF_getprop({P_grid, T_grid}, 'water1');
cases.water1_grid = add_grid_input(out, P_grid, T_grid, [], 'water1');
fprintf('  water1_grid         (%d×%d)\n', numel(P_grid), numel(T_grid));

% ---- NaClaq_5GPa_2024 scatter ------------------------------------------------
P_vec = [100  200  500  900]';
T_vec = [260  270  280  285]';
m_vec = [1.0  2.0  3.0  0.5]';
pts3  = [P_vec, T_vec, m_vec];   % N×3 scatter
out   = SF_getprop(pts3, 'NaClaq_5GPa_2024');
cases.NaClaq_scatter = add_pt_input(out, P_vec, T_vec, m_vec, 'NaClaq_5GPa_2024');
fprintf('  NaClaq_scatter      (%d pts)\n', numel(P_vec));

% ---- NaClaq_5GPa_2024 grid ---------------------------------------------------
P_grid3 = [50 100 300 700];
T_grid3 = [260 270 280 290];
m_grid3 = [0.5 2.0 4.0];
out     = SF_getprop({P_grid3, T_grid3, m_grid3}, 'NaClaq_5GPa_2024');
cases.NaClaq_grid = add_grid_input(out, P_grid3, T_grid3, m_grid3, 'NaClaq_5GPa_2024');
fprintf('  NaClaq_grid         (%d×%d×%d)\n', ...
    numel(P_grid3), numel(T_grid3), numel(m_grid3));

% ---- NaClaq stitched scatter (LP-only, blend, HP-only points) ----------------
P_vec = [200  700  800  2000  5000]';
T_vec = [270  275  280  290   300]';
m_vec = [1.0  2.0  0.5  3.0   1.5]';
pts3  = [P_vec, T_vec, m_vec];
out   = SF_getprop(pts3, 'NaClaq');
cases.NaClaq_stitch_scatter = add_pt_input(out, P_vec, T_vec, m_vec, 'NaClaq');
fprintf('  NaClaq_stitch_scatter (%d pts)\n', numel(P_vec));

% ---- NaClaq stitched grid (P axis spans LP, blend, HP) -----------------------
P_grid4 = [100 400 700 900 1500 3000];
T_grid4 = [260 280 300];
m_grid4 = [0.5 2.0 4.0];
out     = SF_getprop({P_grid4, T_grid4, m_grid4}, 'NaClaq');
cases.NaClaq_stitch_grid = add_grid_input(out, P_grid4, T_grid4, m_grid4, 'NaClaq');
fprintf('  NaClaq_stitch_grid    (%d×%d×%d)\n', ...
    numel(P_grid4), numel(T_grid4), numel(m_grid4));

% ---- water3 (Helmholtz fluid) scatter: vapour, liquid, supercritical ---------
P_vec = [1e-3 0.101325 30 500 2000]';
T_vec = [300  298.15   700 350 450]';
out   = SF_getprop([P_vec, T_vec], 'water3');
cases.water3_scatter = add_pt_input(out, P_vec, T_vec, [], 'water3');
fprintf('  water3_scatter      (%d pts)\n', numel(P_vec));

% ---- water3 grid --------------------------------------------------------------
P_grid = [0.1 10 100 1000];
T_grid = [300 400 500 700];
out    = SF_getprop({P_grid, T_grid}, 'water3');
cases.water3_grid = add_grid_input(out, P_grid, T_grid, [], 'water3');
fprintf('  water3_grid         (%d×%d)\n', numel(P_grid), numel(T_grid));

% ---- (rho,T) input: water3 direct, water1 through SF_rho2P ----------------------
R_vec = [1e-2 997 1100 1250]';
T_vec = [400 298.15 300 350]';
out   = SF_getprop([R_vec, T_vec], 'water3', [], 'input', 'rhoT');
s = add_pt_input(out, R_vec, T_vec, [], 'water3'); s.input_mode = 'rhoT';
cases.water3_rhoT = s;
fprintf('  water3_rhoT         (%d pts)\n', numel(R_vec));
R_grid = [1000 1050 1100];
T_grid = [280 300 330];
out    = SF_getprop({R_grid, T_grid}, 'water1', [], 'input', 'rhoT');
s = add_grid_input(out, R_grid, T_grid, [], 'water1'); s.input_mode = 'rhoT';
cases.water1_rhoT = s;
fprintf('  water1_rhoT         (%d×%d)\n', numel(R_grid), numel(T_grid));

% ---- Save -------------------------------------------------------------------
out_path = fullfile(here, 'reference_getProp.mat');
save(out_path, '-struct', 'cases');
n = numel(fieldnames(cases));
fprintf('\nWrote %d cases to %s\n', n, out_path);
end


% ---- Helper: tag a scatter case with its input PT(m) ------------------------
function s = add_pt_input(out, P, T, m, material)
    s = out;
    s.input_P        = P(:);
    s.input_T        = T(:);
    s.input_type     = 'scatter';
    s.input_material = material;
    if ~isempty(m)
        s.input_m = m(:);
    else
        s.input_m = NaN;
    end
end

% ---- Helper: tag a grid case with its input PT(m) axes ---------------------
function s = add_grid_input(out, P, T, m, material)
    s = out;
    s.input_P        = P(:);
    s.input_T        = T(:);
    s.input_type     = 'grid';
    s.input_material = material;
    if ~isempty(m)
        s.input_m = m(:);
    else
        s.input_m = NaN;
    end
end
