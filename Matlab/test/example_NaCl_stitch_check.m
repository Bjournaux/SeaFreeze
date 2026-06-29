%% example_NaCl_stitch_check.m
% Diagnostic figures for SF_NaCl_stitch: compares the LP spline, HP spline,
% the stitched result, and the older NaClaq spline for all NaCl(aq)
% properties across the blend zone.
%
% For each property a figure is produced with 3×3 subplots:
%   rows = T values [260, 300, 400] K
%   cols = m values [0.5, 2, 5] mol/kg
%
% Each subplot shows:
%   ─── black solid    : stitched (full P range)
%   ─ ─ blue dashed    : LP spline  (P ≤ 1001 MPa)
%   ─ ─ red  dash-dot  : HP spline  (P ≥ 499.9 MPa)
%   ─── magenta dotted : older NaClaq (sp_NaCl_5GPa_500K, full P range)
%   vertical dotted lines at P = 499.9 and 1001 MPa
%
% Smoothness metric (printed to console):
%   max |Δ²f / ΔP²| in the blend zone for each property × (T,m) combo.
%
% Requirements:
%   - Run from Matlab/ or add Matlab/ to the path.
%   - newLBFNaCl/LBF2026.mat and SeaFreeze_Gibbs_VII_NaCl5GPa.mat must be
%     accessible.
%
% Baptiste Journaux — 2026

clear; close all;

%% ---- Paths and data ----------------------------------------------------
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));                  % Matlab/
addpath(fullfile(fileparts(here), 'LocalBasisFunction'));

spLP = sf_load_spline('NaClaq_LP');
spHP = sf_load_spline('NaClaq_HP');

% Older NaClaq spline (currently used by SF_getprop 'NaClaq')
try
    spOLD   = sf_load_spline('NaClaq_5GPa_2024');
    has_old = true;
    fprintf('Older NaClaq: P=[%.1f, %.1f] MPa  T=[%.1f, %.1f] K\n', ...
            spOLD.knots{1}(1), spOLD.knots{1}(end), ...
            spOLD.knots{2}(1), spOLD.knots{2}(end));
catch ME
    warning('example_NaCl_stitch_check:noOldFile', ...
            'Could not load NaClaq spline (%s) — older curve will be skipped.', ME.message);
    spOLD   = [];
    has_old = false;
end

%% ---- Grid parameters ---------------------------------------------------
P_LO  = 499.9;     % MPa — start of blend zone
P_HI  = 1001.0;    % MPa — end of blend zone
dP    = 5;         % MPa — plotting step

P_full = (200 : dP : 2000)';
P_LP   = P_full(P_full <= P_HI);
P_HP   = P_full(P_full >= P_LO);
P_fine = (P_LO : 1 : P_HI)';

% Valid range for older spline (T ≤ 501 K, P ≤ 5000 MPa)
if has_old
    T_old_max = spOLD.knots{2}(end);   % 501 K
    P_old_max = spOLD.knots{1}(end);   % 5000 MPa
    P_old = P_full(P_full <= P_old_max);
else
    T_old_max = -Inf;
    P_old = [];
end

T_vals = [260, 300, 400];   % K
m_vals = [0.5, 2.0, 5.0];  % mol/kg
nT = numel(T_vals);
nm = numel(m_vals);

%% ---- Property groups and metadata --------------------------------------
prop_meta = {
    'G',               'G',                       'J kg^{-1}'
    'S',               'S',                       'J kg^{-1} K^{-1}'
    'U',               'U',                       'J kg^{-1}'
    'H',               'H',                       'J kg^{-1}'
    'A',               'A (Helmholtz)',            'J kg^{-1}'
    'rho',             '\rho',                    'kg m^{-3}'
    'alpha',           '\alpha',                  'K^{-1}'
    'Kt',              'K_T',                     'MPa'
    'Ks',              'K_S',                     'MPa'
    'Kp',              'K''',                     '—'
    'Cp',              'C_P',                     'J kg^{-1} K^{-1}'
    'Cv',              'C_V',                     'J kg^{-1} K^{-1}'
    'vel',             'v_{sound}',               'm s^{-1}'
    'Js',              'J_S (dT/dP|_S)',           'K MPa^{-1}'
    'gamma_Gruneisen', '\gamma_{Gr}',             '—'
    'muw',             '\mu_w',                   'J mol^{-1}'
    'mus',             '\mu_s',                   'J mol^{-1}'
    'phi',             '\phi (osmotic)',           '—'
    'aw',              'a_w (activity)',           '—'
    'Va',              'V_a (app. molar vol)',     'cm^3 mol^{-1}'
    'Cpa',             'C_{pa} (app. mol Cp)',    'J K^{-1} mol^{-1}'
    'Vm',              'V_m (part. mol vol)',      'cm^3 mol^{-1}'
    'Vw',              'V_w (part. mol vol H_2O)','cm^3 mol^{-1}'
    'Cpm',             'C_{pm} (part. mol Cp)',   'J K^{-1} mol^{-1}'
    'Vex',             'V_{ex} (excess vol)',      'cm^3 mol^{-1}'
};

groups = {
    'Thermodynamic potentials',  {'G','S','U','H','A'}
    'Volumetric',                {'rho','alpha','Kt','Ks','Kp'}
    'Thermal & acoustic',        {'Cp','Cv','vel','Js','gamma_Gruneisen'}
    'Chemical (NaCl mixing)',    {'muw','mus','phi','aw'}
    'Apparent / partial molar',  {'Va','Cpa','Vm','Vw','Cpm','Vex'}
};

meta_map = struct();
for r = 1:size(prop_meta,1)
    tag = prop_meta{r,1};
    meta_map.(tag).label = prop_meta{r,2};
    meta_map.(tag).unit  = prop_meta{r,3};
end

%% ---- Pre-evaluate all splines ------------------------------------------
fprintf('Evaluating splines...\n');

res_stitch = cell(nT, nm);
res_LP_raw = cell(nT, nm);
res_HP_raw = cell(nT, nm);
res_OLD    = cell(nT, nm);   % older NaClaq

for iT = 1:nT
    T = T_vals(iT);
    for im = 1:nm
        mval = m_vals(im);
        fprintf('  T=%g K, m=%g mol/kg\n', T, mval);

        res_stitch{iT,im} = SF_NaCl_stitch(spLP, spHP, {P_full, T, mval});
        res_LP_raw{iT,im} = fnGval(spLP, {P_LP, T, mval});
        res_HP_raw{iT,im} = fnGval(spHP, {P_HP, T, mval});

        if has_old && T <= T_old_max
            res_OLD{iT,im} = fnGval(spOLD, {P_old, T, mval});
        end
    end
end

%% ---- Smoothness check --------------------------------------------------
fprintf('\n=== Smoothness check (max |d²f/dP²| in blend zone [%.1f, %.1f] MPa) ===\n', ...
        P_LO, P_HI);
fprintf('%-22s', 'Property');
for iT = 1:nT
    for im = 1:nm
        fprintf('  T=%3g m=%3.1f', T_vals(iT), m_vals(im));
    end
end
fprintf('\n%s\n', repmat('-', 1, 22 + (nT*nm)*14));

all_tags = prop_meta(:,1);
for ip = 1:numel(all_tags)
    tag = all_tags{ip};
    fprintf('%-22s', tag);
    for iT = 1:nT
        T = T_vals(iT);
        for im = 1:nm
            mval = m_vals(im);
            try
                r_fine = SF_NaCl_stitch(spLP, spHP, {P_fine, T, mval}, tag);
                vals = squeeze(r_fine.(tag));
                d2   = diff(diff(vals));
                fprintf('  %12.3g', max(abs(d2)));
            catch
                fprintf('  %12s', 'n/a');
            end
        end
    end
    fprintf('\n');
end
fprintf('\n');

%% ---- Plotting ----------------------------------------------------------
clr_stitch = [0    0    0   ];   % black solid
clr_LP     = [0    0.4  0.8 ];   % blue dashed
clr_HP     = [0.8  0.1  0.1 ];   % red dash-dot
clr_old    = [0.6  0    0.8 ];   % magenta/purple dotted

lw_stitch = 2.0;
lw_raw    = 1.2;
lw_old    = 1.4;

for ig = 1:size(groups, 1)
    grp_title = groups{ig, 1};
    grp_tags  = groups{ig, 2};

    for ip = 1:numel(grp_tags)
        tag = grp_tags{ip};
        if ~isfield(meta_map, tag), continue; end

        fig = figure('Name', sprintf('%s — %s', grp_title, tag), ...
                     'NumberTitle', 'off', ...
                     'Units', 'normalized', ...
                     'Position', [0.05 0.05 0.88 0.85]);

        ha = gobjects(nT, nm);
        for iT = 1:nT
            for im = 1:nm
                ha(iT,im) = subplot(nT, nm, (iT-1)*nm + im);
            end
        end

        for iT = 1:nT
            T    = T_vals(iT);
            for im = 1:nm
                mval = m_vals(im);
                ax   = ha(iT, im);

                % --- extract data (squeeze to 1-D) ----------------------
                try
                    y_s  = squeeze(res_stitch{iT,im}.(tag));
                    y_lp = squeeze(res_LP_raw{iT,im}.(tag));
                    y_hp = squeeze(res_HP_raw{iT,im}.(tag));
                catch
                    text(ax, 0.5, 0.5, 'N/A', ...
                         'Units','normalized','HorizontalAlignment','center');
                    title(ax, sprintf('T=%g K, m=%g mol/kg', T, mval));
                    continue
                end

                % Older NaClaq: available only within its T range
                has_old_here = has_old && T <= T_old_max && ...
                               ~isempty(res_OLD{iT,im}) && ...
                               isfield(res_OLD{iT,im}, tag);
                if has_old_here
                    y_old = squeeze(res_OLD{iT,im}.(tag));
                else
                    y_old = [];
                end

                axes(ax); %#ok<LAXES>
                hold on;

                % --- plot curves ----------------------------------------
                h_s  = plot(P_full, y_s,  '-',  'Color', clr_stitch, 'LineWidth', lw_stitch);
                h_lp = plot(P_LP,   y_lp, '--', 'Color', clr_LP,     'LineWidth', lw_raw);
                h_hp = plot(P_HP,   y_hp, '-.', 'Color', clr_HP,     'LineWidth', lw_raw);

                leg_h = [h_s, h_lp, h_hp];
                leg_l = {'Stitched', 'LP spline', 'HP spline'};

                if has_old_here
                    h_old = plot(P_old, y_old, ':', ...
                                 'Color', clr_old, 'LineWidth', lw_old);
                    leg_h(end+1) = h_old; %#ok<AGROW>
                    leg_l{end+1} = 'older NaClaq'; %#ok<AGROW>
                end

                % --- blend-zone boundary lines --------------------------
                all_y = [y_s; y_lp; y_hp];
                if has_old_here, all_y = [all_y; y_old]; end %#ok<AGROW>
                yl_pad = [min(all_y), max(all_y)];
                if diff(yl_pad) == 0, yl_pad = yl_pad + [-1 1]; end
                plot([P_LO P_LO], yl_pad, 'k:', 'LineWidth', 0.8, 'HandleVisibility','off');
                plot([P_HI P_HI], yl_pad, 'k:', 'LineWidth', 0.8, 'HandleVisibility','off');

                box on;
                xlabel('P (MPa)');
                ylabel(sprintf('%s (%s)', meta_map.(tag).label, meta_map.(tag).unit), ...
                       'Interpreter', 'tex');
                title(ax, sprintf('T=%g K,  m=%g mol/kg', T, mval));

                if iT == 1 && im == 1
                    legend(leg_h, leg_l, 'Location', 'best', 'FontSize', 7);
                end
            end
        end

        sgtitle(sprintf('%s:  %s', grp_title, meta_map.(tag).label), ...
                'Interpreter', 'tex', 'FontSize', 13, 'FontWeight', 'bold');
    end
end

fprintf('Done — %d figures created.\n', numel(findobj('Type','figure')));
