function [fig, pm] = SF_WPD_PT(varargin)
% SF_WPD_PT  Full H2O phase diagram in (P,T): vapour, liquid, supercritical
% fluid, critical point and the ices, by Gibbs-energy minimisation.
% Baptiste Journaux - 2026
%
% The fluid is the Helmholtz material water3 (vapour + liquid + supercritical),
% the ices are the SeaFreeze Gibbs splines.  At every grid point the stable
% phase is the one with the lowest specific Gibbs energy (sf_phase_map);
% phase boundaries are the G_i = G_j contours between neighbouring stable
% phases, the saturation curve comes from SF_coexistence.
%
% Usage:
%   SF_WPD_PT()                                   % defaults below
%   SF_WPD_PT('P', [1e-8 1e5], 'T', [150 1800])   % MPa, K (log P axis)
%   SF_WPD_PT('ices', {'Ih','II','III','V','VI','VII_X_French'})
%   SF_WPD_PT('ax', gca, 'nP', 300, 'nT', 250)
%   [fig, pm] = SF_WPD_PT(...)                    % pm: the phase map (sf_phase_map)
%
% Name-value parameters:
%   'ax'     axes to draw into (default: new figure)
%   'P'      [Pmin Pmax] MPa, log-spaced            (default [1e-8 1e5])
%   'T'      [Tmin Tmax] K                          (default [150 1800])
%   'nP','nT' grid size                             (default 400 x 320)
%   'fluid'  Helmholtz fluid                        (default 'water3')
%   'ices'   ice phases                             (default Ih, II, III, V, VI;
%            ice VII/X is left out until an updated model is available)
%
% Where no phase is available (the ice VII/X field when VII is not
% included, below the fluid's validity mask) the diagram is grey and marked
% "not modelled".  See sf_phase_map for the validity masks.
%
% See also: SF_WPD_rhoT, SF_WPD, SF_PhaseLines, SF_coexistence.

p = inputParser;
addParameter(p, 'ax', []);
addParameter(p, 'P', [1e-8 1e5]);
addParameter(p, 'T', [150 1800]);
addParameter(p, 'nP', 400);
addParameter(p, 'nT', 320);
addParameter(p, 'fluid', 'water3');
addParameter(p, 'ices', {'Ih','II','III','V','VI'});
parse(p, varargin{:});
o = p.Results;
warning('SeaFreeze:longRuntime', ['SF_WPD_PT: computing the stable phase at %d (P,T) states ' ...
        '(Gibbs minimisation over the Helmholtz fluid and the ices). This typically takes ' ...
        '0.5-2 min depending on your machine; reduce ''nP''/''nT'' for a faster preview. ' ...
        'Silence with warning(''off'',''SeaFreeze:longRuntime'').'], o.nP * o.nT);

Pg = logspace(log10(o.P(1)), log10(o.P(2)), o.nP)';
Tg = linspace(o.T(1), o.T(2), o.nT);
pm = sf_phase_map(Pg, Tg, 'fluid', o.fluid, 'ices', o.ices);
n = numel(pm.names);

if isempty(o.ax)
    fig = figure('Position', [100 100 1000 680]); ax = axes(fig);
else
    ax = o.ax; fig = ancestor(ax, 'figure');
end
hold(ax, 'on');
[cols, cmiss] = sf_wpd_colors(n);
% fields: an RGB image on the log-P grid (NaN-free, works in MATLAB and Octave)
img = sf_wpd_image(pm.stable, cols, cmiss);           % nP x nT x 3
image(ax, 'XData', log10(Pg), 'YData', Tg, 'CData', permute(img, [2 1 3]));
set(ax, 'YDir', 'normal');
% boundaries between neighbouring stable phases
for i = 1:n
    for j = i+1:n
        if ~any(pm.stable(:) == i) || ~any(pm.stable(:) == j), continue; end
        Z = reshape(pm.G(i, :, :) - pm.G(j, :, :), numel(Pg), numel(Tg));
        Z(~(pm.stable == i | pm.stable == j)) = NaN;
        segs = sf_contour_segments(log10(Pg), Tg, Z.');
        for k = 1:numel(segs)
            plot(ax, segs{k}(:, 1), segs{k}(:, 2), 'k-', 'LineWidth', 1.1);
        end
    end
end
% vapour-liquid saturation curve (stable branch) and critical point
Tc = 647.096; Pc = 22.064;
Ts = unique([linspace(max(273.16, o.T(1)), Tc - 1, 80), Tc - logspace(log10(0.05), log10(30), 40)]);
sat = SF_coexistence('saturation', Ts, 'fluid', o.fluid);
plot(ax, log10(sat.P), sat.T, 'k-', 'LineWidth', 1.1);
plot(ax, log10(Pc), Tc, 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 1.4);
text(ax, log10(Pc) + 0.1, Tc - 25, 'critical point', 'FontSize', 9);
% triple points
tp = sf_triple_points(pm, o.fluid);
for k = 1:numel(tp)
    plot(ax, log10(tp(k).P), tp(k).T, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
end
% labels
sf_wpd_labels(ax, log10(Pg), Tg, pm, 'PT');
xlim(ax, log10(o.P)); ylim(ax, o.T);
xt = ceil(log10(o.P(1))):2:floor(log10(o.P(2)));
set(ax, 'XTick', xt, 'XTickLabel', arrayfun(@(e) sprintf('10^{%d}', e), xt, 'UniformOutput', false));
xlabel(ax, 'Pressure (MPa)'); ylabel(ax, 'Temperature (K)');
title(ax, sprintf('H_2O phase diagram — fluid: %s, ices: %s', o.fluid, strjoin(strrep(o.ices, 'VII_X_French', 'VII/X'), ', ')), ...
      'Interpreter', 'tex');
box(ax, 'on');
if nargout == 0, clear fig; end
end
