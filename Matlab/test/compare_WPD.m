% compare_WPD.m — run from Matlab/ directory
% Visual comparison of old SF_WPD_v1 vs new SF_WPD, plus timing.
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..');

% --- Side-by-side comparison ---
fig1 = figure('Position', [100 100 1400 600], 'Visible', 'off');
subplot(1,2,1);
SF_WPD_v1();
title('SF\_WPD v1 (old - static WPD.mat)', 'FontSize', 12);

ax2 = subplot(1,2,2);
SF_WPD('ax', ax2, 'labels', true);
title(ax2, 'SF\_WPD v2 (new - dynamic)', 'FontSize', 12);

saveas(fig1, fullfile(here, 'compare_WPD_v1_v2.png'));
fprintf('Saved test/compare_WPD_v1_v2.png\n');

% --- NaClaq overlay ---
fig2 = SF_WPD('solute', 'NaCl', 'm', [0.5 1 2 4], 'labels', true);
set(fig2, 'Visible', 'off');
saveas(fig2, fullfile(here, 'compare_WPD_NaCl.png'));
fprintf('Saved test/compare_WPD_NaCl.png\n');

% --- Timing comparison ---
fprintf('\n--- Timing ---\n');
clear functions;
tic; SF_WPD_v1(); t_old = toc; close;
fprintf('Old (static):     %.3f s\n', t_old);

clear functions;
tic; SF_WPD(); t_cold = toc; close;
fprintf('New (cold):       %.3f s\n', t_cold);

tic; SF_WPD(); t_warm = toc; close;
fprintf('New (warm cache): %.3f s\n', t_warm);

close all;
