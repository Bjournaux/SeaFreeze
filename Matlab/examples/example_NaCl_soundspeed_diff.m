% example_NaCl_soundspeed_diff
%
% Plot the relative sound-speed difference (in %) between aqueous NaCl
% solutions and pure water as a function of pressure, at fixed temperature,
% for three molalities.  Six isotherms are overlaid: T = 30–80 °C.
%
% Difference plotted (PPM):
%     dvel_pct(P, T, m) = [vel_NaClaq(P, T, m) - vel_water1(P, T)]
%                          / vel_water1(P, T)  ×  100
%
% Uses water1 (Bollengier et al. 2019) as the pure-water reference.
%
% Baptiste Journaux - 2026

%% --- Parameters ---------------------------------------------------------
P  = 0.1:1:700;                % Pressure (MPa)
Ts = [303.15, 313.15, 323.15, 333.15, 343.15, 353.15]; % Temperatures (K) = 30–80 °C
ms = [1.0, 3.0, 6.0];          % Molalities (mol/kg)

%% --- Add SeaFreeze to path ----------------------------------------------
% LocalBasisFunction is added first so the (top-level) Matlab/ folder ends
% up at the front of the MATLAB path. The modern fnGval/sp_val live there
% and need to shadow the older copies inside LocalBasisFunction/.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'LocalBasisFunction'));
addpath(fullfile(here, '..'));

%% --- Compute -----------------------------------------------------------
diffs = zeros(numel(P), numel(ms), numel(Ts));
for it = 1:numel(Ts)
    T = Ts(it);
    out_w  = SF_getprop({P, T}, 'water1', 'vel');
    vel_w  = out_w.vel(:);                     % (N_P x 1)
    for im = 1:numel(ms)
        out_n = SF_getprop({P, T, ms(im)}, 'NaClaq', 'vel');
        % NaClaq returns (N_P x 1 x 2) for a scalar m because the internal
        % baseline slice at m = cutoff is not stripped when only base
        % properties are requested; take the last slice = user's m.
        vel_n = out_n.vel(:, :, end);          % (N_P x 1)
        diffs(:, im, it) = (vel_n(:) - vel_w) ./ vel_w * 100;   % %
    end
end

%% --- Plot --------------------------------------------------------------
% Encoding: colour = temperature (cool colormap), line style = molality.
% A compact legend shows one swatch per T and one per m rather than
% listing all numel(Ts)*numel(ms) curves individually.

figure('Position', [50 50 960 620], 'Color', 'w');
hold on; grid on;

cmap   = cool(numel(Ts));           % one colour per temperature
styles = {'-', '--', ':'};          % one line style per molality (up to 3)

for it = 1:numel(Ts)
    for im = 1:numel(ms)
        plot(P, diffs(:, im, it), styles{im}, ...
             'Color', cmap(it, :), 'LineWidth', 1.5, ...
             'HandleVisibility', 'off');
    end
end

yline(0, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
xlabel('Pressure (MPa)');
ylabel('Difference in Sound Speeds  (%)');
title('Sound-speed shift due to NaCl');

% --- compact legend: colour key for T, style key for m ----------------
ax = gca;
h_T = gobjects(numel(Ts), 1);
for it = 1:numel(Ts)
    h_T(it) = plot(nan, nan, '-', 'Color', cmap(it,:), 'LineWidth', 2.0, ...
        'DisplayName', sprintf('T = %g °C', Ts(it) - 273.15));
end
h_m = gobjects(numel(ms), 1);
for im = 1:numel(ms)
    h_m(im) = plot(nan, nan, styles{im}, 'Color', [0.15 0.15 0.15], 'LineWidth', 1.5, ...
        'DisplayName', sprintf('m = %g mol/kg', ms(im)));
end
legend([h_T; h_m], 'Location', 'best', 'FontSize', 9);

xlim([min(P), max(P)]);
hold off;
