function [fig, pm] = SF_WPD_rhoT(varargin)
% SF_WPD_rhoT  Full H2O phase diagram in (rho,T): vapour, liquid,
% supercritical fluid and the ices, with the two-phase (coexistence)
% regions — vapour-liquid dome, sublimation, melting and ice-ice.
% Baptiste Journaux - 2026
%
% The (P,T) phase map of sf_phase_map is re-drawn in density: along every
% isotherm the density of the stable phase increases with P; where it jumps
% (between two phases, or across the vapour-liquid saturation) the density
% gap is a two-phase region, drawn grey and bounded by the coexisting
% densities.  The vapour-liquid dome is bounded exactly by the saturation
% densities from SF_coexistence.
%
% Usage:
%   SF_WPD_rhoT()                                           % log-density view
%   SF_WPD_rhoT('rho', [850 1700], 'T', [150 500], 'P', [1e-10 1e4], 'xscale', 'linear')
%   [fig, pm] = SF_WPD_rhoT(...)
%
% Name-value parameters:
%   'ax'      axes to draw into (default: new figure)
%   'rho'     [rhomin rhomax] kg/m^3              (default [1e-7 4e3])
%   'T'       [Tmin Tmax] K                       (default [150 1800])
%   'P'       [Pmin Pmax] MPa of the underlying P-T map (default [1e-10 1e5])
%   'nP','nT' P-T grid size                       (default 900 x 320)
%   'nrho'    density pixels                      (default 700)
%   'xscale'  'log' (default) or 'linear'
%   'fluid', 'ices'  as SF_WPD_PT
%
% See also: SF_WPD_PT, sf_phase_map, SF_coexistence.

p = inputParser;
addParameter(p, 'ax', []);
addParameter(p, 'rho', [1e-7 4e3]);
addParameter(p, 'T', [150 1800]);
addParameter(p, 'P', [1e-10 1e5]);
addParameter(p, 'nP', 900);
addParameter(p, 'nT', 320);
addParameter(p, 'nrho', 700);
addParameter(p, 'xscale', 'log');
addParameter(p, 'fluid', 'water_Brown2026');
addParameter(p, 'ices', {'Ih','II','III','V','VI'});
parse(p, varargin{:});
o = p.Results;
o.fluid = sf_material_name(char(o.fluid));   % renamed materials
warning('SeaFreeze:longRuntime', ['SF_WPD_rhoT: computing the stable phase at %d (P,T) states ' ...
        '(Gibbs minimisation over the Helmholtz fluid and the ices). This typically takes ' ...
        '0.5-2 min depending on your machine; reduce ''nP''/''nT'' for a faster preview. ' ...
        'Silence with warning(''off'',''SeaFreeze:longRuntime'').'], o.nP * o.nT);
islog = strcmpi(o.xscale, 'log');

Pg = logspace(log10(o.P(1)), log10(o.P(2)), o.nP)';
Tg = linspace(o.T(1), o.T(2), o.nT);
pm = sf_phase_map(Pg, Tg, 'fluid', o.fluid, 'ices', o.ices);
n = numel(pm.names); TWO = n + 1;
if islog
    rq = logspace(log10(o.rho(1)), log10(o.rho(2)), o.nrho)'; X = log10(rq);
else
    rq = linspace(o.rho(1), o.rho(2), o.nrho)'; X = rq;
end
fx = @(r) r; if islog, fx = @(r) log10(r); end

% saturation densities (solved on ~100 T clustered toward Tc, interpolated)
Tc = 647.096; rhoc = 322;
Ts = unique([linspace(273.16, Tc - 1, 80), Tc - logspace(log10(0.05), log10(30), 40)]);
sat = SF_coexistence('saturation', Ts, 'fluid', o.fluid);
okS = isfinite(sat.P);

img = zeros(o.nrho, o.nT);
acc = zeros(n, n, 3);                      % per pair: sum X, sum T, count
for j = 1:o.nT
    st = pm.stable(:, j); rs = pm.rho_stable(:, j);
    ok = st > 0 & isfinite(rs);
    if nnz(ok) < 2, continue; end
    st = st(ok); rs = cummax(rs(ok));
    k = sum(bsxfun(@lt, rs(:).', rq), 2);           % rs(k) < rho <= rs(k+1)
    inside = k > 0 & k < numel(rs);
    kk = min(max(k, 1), numel(rs) - 1);
    a = st(kk); b = st(kk + 1);
    single = inside & a == b;
    dual = inside & a ~= b;
    img(single, j) = a(single);
    img(dual, j) = TWO;
    lo = min(a(dual), b(dual)); hi = max(a(dual), b(dual));
    idd = find(dual);
    for q = 1:numel(idd)
        acc(lo(q), hi(q), 1) = acc(lo(q), hi(q), 1) + X(idd(q));
        acc(lo(q), hi(q), 2) = acc(lo(q), hi(q), 2) + Tg(j);
        acc(lo(q), hi(q), 3) = acc(lo(q), hi(q), 3) + 1;
    end
    % vapour-liquid dome
    if Tg(j) >= 273.16 && Tg(j) < Tc && any(okS)
        rv = exp(interp1(Ts(okS), log(sat.rho_B(okS)), Tg(j)));
        rl = interp1(Ts(okS), sat.rho_A(okS), Tg(j));
        dome = rq > rv & rq < rl & img(:, j) == 1;
        img(dome, j) = TWO;
        acc(1, 1, :) = acc(1, 1, :) + reshape([sum(X(dome)) Tg(j) * nnz(dome) nnz(dome)], 1, 1, 3);
    end
end

if isempty(o.ax)
    fig = figure('Position', [100 100 1000 680]); ax = axes(fig);
else
    ax = o.ax; fig = ancestor(ax, 'figure');
end
hold(ax, 'on');
[cols, ~, ctwo] = sf_wpd_colors(n);
im = sf_wpd_image(img, cols, [1 1 1], ctwo);
image(ax, 'XData', X, 'YData', Tg, 'CData', permute(im, [2 1 3]));
set(ax, 'YDir', 'normal');

% coexisting densities along every P-T boundary
sp = sf_load_spline(o.fluid);
for i = 1:n
    for jj = i+1:n
        if ~any(pm.stable(:) == i) || ~any(pm.stable(:) == jj), continue; end
        Z = reshape(pm.G(i, :, :) - pm.G(jj, :, :), numel(Pg), numel(Tg));
        Z(~(pm.stable == i | pm.stable == jj)) = NaN;
        segs = sf_contour_segments(log10(Pg), Tg, Z.');
        for s = 1:numel(segs)
            Pl = 10.^segs{s}(:, 1); Tl = segs{s}(:, 2);
            for ph = [i jj]
                r = rho_on(pm.names{ph}, ph == 1, sp, Pl, Tl);
                plot(ax, fx(r), Tl, 'k-', 'LineWidth', 0.8);
            end
        end
    end
end
plot(ax, fx(sat.rho_A(okS)), Ts(okS), 'k-', 'LineWidth', 1.1);
plot(ax, fx(sat.rho_B(okS)), Ts(okS), 'k-', 'LineWidth', 1.1);
plot(ax, fx(rhoc), Tc, 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'LineWidth', 1.4);
text(ax, fx(rhoc), Tc + 30, 'critical point', 'FontSize', 9, 'HorizontalAlignment', 'right');
% three-phase (triple-point) tie lines: the separations between two-phase regions
tp = sf_triple_points(pm, o.fluid);
for k = 1:numel(tp)
    r = tp(k).rho;
    plot(ax, fx([min(r) max(r)]), [tp(k).T tp(k).T], 'k-', 'LineWidth', 1.0);
    plot(ax, fx(r), tp(k).T * [1 1 1], 'k|', 'MarkerSize', 8, 'LineWidth', 1.2);
end

% labels: single-phase fields and two-phase regions
pm.xrhoc = fx(rhoc);
sf_wpd_labels(ax, X, Tg, pm, 'rhoT', img);
lab = strrep(pm.names, 'VII_X_French', 'VII/X');
for i = 1:n
    for jj = i:n
        c = acc(i, jj, 3);
        if c < 60, continue; end
        xm = acc(i, jj, 1) / c; tm = acc(i, jj, 2) / c;
        if i == 1 && jj == 1
            s = 'L + V';
        elseif i == 1
            rm = xm; if islog, rm = 10^xm; end
            if rm < 100, s = [lab{jj} ' + V']; else, s = [lab{jj} ' + L']; end
        else
            s = [lab{i} ' + ' lab{jj}];
        end
        text(ax, xm, tm, s, 'FontSize', 8.5, 'FontAngle', 'italic', 'Color', [0.32 0.32 0.31], ...
             'HorizontalAlignment', 'center');
    end
end
xlim(ax, [X(1) X(end)]); ylim(ax, o.T);
if islog
    xt = ceil(log10(o.rho(1))):2:floor(log10(o.rho(2)));
    set(ax, 'XTick', xt, 'XTickLabel', arrayfun(@(e) sprintf('10^{%d}', e), xt, 'UniformOutput', false));
end
xlabel(ax, 'Density (kg/m^3)'); ylabel(ax, 'Temperature (K)');
title(ax, sprintf('H_2O phase diagram in density — fluid: %s; grey: two-phase regions', o.fluid));
box(ax, 'on');
if nargout == 0, clear fig; end
end


function r = rho_on(mat, isfluid, sp, P, T)
% density of one phase at scattered (P,T) points on a boundary
o = SF_getprop([P(:) T(:)], mat, 'rho');
r = o.rho(:);
if isfluid
    lo = ~isfinite(r);
    if any(lo), [~, r(lo)] = sf_ideal_gas(sp, P(lo), T(lo)); end
end
end
