% compare_WPD_overlap.m — overlay old (black) and new (dashed red)
here = fileparts(mfilename('fullpath'));

fig = figure('Position', [100 100 900 650], 'Visible', 'off');
ax = axes(fig); hold(ax, 'on');

% Old v1 — black solid/dashed (plots into current axes)
SF_WPD_v1();
% Hide v1 plot lines from legend
ch = get(ax, 'Children');
set(ch, 'HandleVisibility', 'off');

% New v2 — dashed red, overlaid on same axes (stable only)
PAIRS = {'Ih','water1'; 'Ih','II'; 'Ih','III'; ...
         'II','III'; 'II','V'; 'II','VI'; 'II','water1'; ...
         'III','V'; 'III','water1'; ...
         'V','water1'; 'V','VI'; 'VI','water1'};
for k = 1:size(PAIRS, 1)
    try
        r = SF_PhaseLines(PAIRS{k,1}, PAIRS{k,2}, 'segment', 'all');
    catch
        continue
    end
    if any(r.stable)
        plot(ax, r.P(r.stable), r.T(r.stable), 'r--', 'LineWidth', 1.5, ...
             'HandleVisibility', 'off');
    end
end

% Dummy lines for legend
plot(ax, NaN, NaN, 'k-', 'LineWidth', 1.2, 'DisplayName', 'v1 (old, static)');
plot(ax, NaN, NaN, 'r--', 'LineWidth', 1.5, 'DisplayName', 'v2 (new, dynamic)');
legend(ax, 'show', 'Location', 'southeast', 'FontSize', 11);

title(ax, 'SF\_WPD overlap: v1 (black) vs v2 (dashed red)', 'FontSize', 13);
xlim(ax, [0 2300]); ylim(ax, [0 400]);

saveas(fig, fullfile(here, 'compare_WPD_overlap.png'));
fprintf('Saved test/compare_WPD_overlap.png\n');
close(fig);
