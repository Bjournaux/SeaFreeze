function compare_SF_PhaseLines(out_dir)
% compare_SF_PhaseLines  Generate v1-vs-new comparison figures.
%
% Produces a 4x3 multi-panel PNG showing each pure-phase pair (v1 in black,
% new full curve in red, new stable highlighted, triple points marked) plus
% a NaClaq-melt panel showing the new code's NaCl freezing curves at several
% molalities.
%
% Usage:
%   compare_SF_PhaseLines              % saves to current dir
%   compare_SF_PhaseLines('/tmp')      % saves to /tmp

if nargin < 1, out_dir = pwd; end

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'LocalBasisFunction'));
addpath(fileparts(here));

warning('off','SeaFreeze:nonPositiveT');
warning('off','SeaFreeze:negativePressure');

% =========================================================================
% Figure 1: pure pairs — 11 panels in a 4x3 grid
% =========================================================================
pure_pairs = {
    {'Ih','water1'}, {'Ih','II'},     {'Ih','III'}, ...
    {'II','III'},    {'II','V'},      {'II','VI'}, ...
    {'III','V'},     {'III','water1'},{'V','water1'}, ...
    {'V','VI'},      {'VI','water1'}};

f1 = figure('Position', [50 50 1400 1100], 'Color', 'w');
sgtitle(f1, 'SF\_PhaseLines: v1 vs new (pure phases)', 'FontWeight','bold');

for k = 1:length(pure_pairs)
    p = pure_pairs{k};
    subplot(4, 3, k); hold on; grid on;

    % v1 reference (full curve, black)
    v1 = SF_PhaseLines_v1(p{1}, p{2});
    if ~isempty(v1)
        plot(v1(:,1), v1(:,2), '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 2.0, ...
             'DisplayName','v1 (reference)');
    end

    % new code: full curve (light red) plus stable portion (solid red)
    new_full = SF_PhaseLines(p{1}, p{2}, 'segment','all');
    new_stab = SF_PhaseLines(p{1}, p{2}, 'segment','stable');
    if ~isempty(new_full.P)
        plot(new_full.P, new_full.T, ':r', 'LineWidth', 1.0, 'DisplayName','new (full)');
    end
    if ~isempty(new_stab.P)
        plot(new_stab.P, new_stab.T, '-r', 'LineWidth', 1.4, 'DisplayName','new (stable)');
    end

    % triple points
    if ~isempty(new_full.triple_points)
        plot(new_full.triple_points(:,2), new_full.triple_points(:,1), ...
             'bo', 'MarkerFaceColor','b', 'MarkerSize', 5, 'HandleVisibility','off');
    end

    title(sprintf('%s — %s', p{1}, p{2}), 'Interpreter','none');
    xlabel('P (MPa)'); ylabel('T (K)');
    if k == 1, legend('Location','best','FontSize', 7); end
end

png1 = fullfile(out_dir, 'SF_PhaseLines_v1_vs_new_pure.png');
exportgraphics(f1, png1, 'Resolution', 150);
fprintf('Wrote %s\n', png1);
close(f1);

% =========================================================================
% Figure 2: NaClaq melting curves at multiple molalities
% =========================================================================
ice_phases = {'Ih','III','V','VI'};   % phases with substantial NaClaq curves
molalities = [0.1 0.5 1.0 2.0 4.0];
colors     = lines(length(molalities));

f2 = figure('Position',[50 50 1200 900],'Color','w');
sgtitle(f2, 'SF\_PhaseLines: NaClaq melting curves (new code)', 'FontWeight','bold');

for k = 1:length(ice_phases)
    ice = ice_phases{k};
    subplot(2, 2, k); hold on; grid on;

    % pure-water reference (black, from new code)
    pure = SF_PhaseLines(ice, 'water1', 'segment','all');
    if ~isempty(pure.P)
        plot(pure.P, pure.T, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.8, ...
             'DisplayName','pure water');
    end

    for j = 1:length(molalities)
        m = molalities(j);
        try
            o = SF_PhaseLines(ice, 'NaClaq', 'm', m, 'segment','all');
            if ~isempty(o.P)
                plot(o.P, o.T, '-', 'Color', colors(j,:), 'LineWidth', 1.2, ...
                     'DisplayName', sprintf('m=%.2g mol/kg', m));
            end
        catch err
            fprintf('  %s @ m=%g: %s\n', ice, m, err.message);
        end
    end

    title(sprintf('%s ↔ NaClaq', ice));
    xlabel('P (MPa)'); ylabel('T (K)');
    legend('Location','best','FontSize',8);
end

png2 = fullfile(out_dir, 'SF_PhaseLines_NaCl_melting_curves.png');
exportgraphics(f2, png2, 'Resolution', 150);
fprintf('Wrote %s\n', png2);
close(f2);

% =========================================================================
% Figure 3: full water phase diagram overlay (v1 in grey, new in red)
% =========================================================================
f3 = figure('Position',[50 50 900 700],'Color','w'); hold on; grid on;
title('Water phase diagram: v1 (black) vs new (red)','FontWeight','bold');

for k = 1:length(pure_pairs)
    p = pure_pairs{k};
    v1 = SF_PhaseLines_v1(p{1}, p{2});
    if ~isempty(v1)
        plot(v1(:,1), v1(:,2), '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5, 'HandleVisibility','off');
    end
    n = SF_PhaseLines(p{1}, p{2}, 'segment','stable');
    if ~isempty(n.P)
        plot(n.P, n.T, '-r', 'LineWidth', 1.2, 'HandleVisibility','off');
    end
end

xlabel('P (MPa)'); ylabel('T (K)');
xlim([0 2300]); ylim([180 380]);
text(50, 200, 'Ih', 'FontSize', 12);
text(400, 220, 'II', 'FontSize', 12);
text(280, 250, 'III', 'FontSize', 12);
text(525, 250, 'V', 'FontSize', 12);
text(1000, 250, 'VI', 'FontSize', 12);
text(800, 320, 'Liquid', 'FontSize', 12);

png3 = fullfile(out_dir, 'SF_PhaseLines_full_diagram.png');
exportgraphics(f3, png3, 'Resolution', 150);
fprintf('Wrote %s\n', png3);
close(f3);

fprintf('\nAll comparison figures written to %s\n', out_dir);
end
