function out = SF_PhaseLines(matA, matB, varargin)
% SF_PhaseLines  Equilibrium line between two phases.
%
% Computes the (P, T) curve where the chemical potentials of two phases
% balance. The default sampling grid is the intersection of the two
% phases' spline domains (read from sp.knots, via SF_phase_range).
%
% Usage:
%   out = SF_PhaseLines(matA, matB)                    % full curve, default grid
%   out = SF_PhaseLines(matA, matB, 'm', 1.0)          % NaClaq pair: molality (mol/kg)
%   out = SF_PhaseLines(matA, matB, 'P', 0.1:1:300)    % override pressure grid (MPa)
%   out = SF_PhaseLines(matA, matB, 'T', 240:0.5:280)  % override temperature grid (K)
%   out = SF_PhaseLines(matA, matB, 'segment','stable')% only stable portion
%   out = SF_PhaseLines(matA, matB, 'segment','meta')  % only metastable extension
%   out = SF_PhaseLines(matA, matB, 'plot', true)      % render the curve
%
% Returns a struct with:
%   out.matA, out.matB : phase names
%   out.m              : molality (NaClaq pairs only)
%   out.P  (Nx1, MPa)  : pressure along the equilibrium curve
%   out.T  (Nx1, K)    : temperature along the equilibrium curve
%   out.stable (Nx1)   : logical, true on the thermodynamically stable portion
%   out.segment        : 'all' | 'stable' | 'meta' (the segment selection applied)
%   out.fig            : figure handle (only present when 'plot', true)
%
% For ice ↔ NaClaq pairs the equilibrium condition is the chemical-potential
% balance for water: G_ice * MW_H2O = muw_NaClaq.
%
% Stable / metastable classification uses fixed per-pair triple-point bounds
% (see PAIRS_TABLE below); this matches v1 SF_PhaseLines (now SF_PhaseLines_v1)
% by design — no auto-comparison against other phases.
%
% Baptiste Journaux 2026 (rewrite of Clinton & Journaux 2020).

% ------------ Parse inputs --------------------------------------------------
p = inputParser;
addRequired(p, 'matA', @(s) ischar(s) || (isstring(s) && isscalar(s)));
addRequired(p, 'matB', @(s) ischar(s) || (isstring(s) && isscalar(s)));
addParameter(p, 'm',       NaN,    @(x) isnumeric(x) && (isscalar(x) || isvector(x)));
addParameter(p, 'P',       [],     @(x) isnumeric(x) && isvector(x));
addParameter(p, 'T',       [],     @(x) isnumeric(x) && isvector(x));
addParameter(p, 'segment', 'all',  @(s) ismember(lower(char(s)), {'all','stable','meta'}));
addParameter(p, 'plot',    false,  @valid_plot_arg);
parse(p, matA, matB, varargin{:});
matA = char(matA);
matB = char(matB);
m_val   = p.Results.m(:)';   % row vector (scalar stays scalar)
P_user  = p.Results.P;
T_user  = p.Results.T;
segment = lower(char(p.Results.segment));
% 'plot' may be a logical/scalar (true => new figure, false => no plot),
% or an existing figure handle to overlay onto.
plot_arg   = p.Results.plot;
if isgraphics(plot_arg, 'figure')
    do_plot    = true;
    target_fig = plot_arg;
else
    do_plot    = logical(plot_arg);
    target_fig = [];
end

% ------------ Look up the pair (handles symmetric ordering) -----------------
[pair, swapped] = lookup_pair(matA, matB);
if swapped, [matA, matB] = deal(matB, matA); end

% ------------ Validate molality if NaClaq is involved -----------------------
nacl_involved = any(strcmp({matA, matB}, 'NaClaq'));
if nacl_involved
    if any(isnan(m_val))
        error('SeaFreeze:badInput', ...
            'Pair (%s, %s) involves NaClaq; pass molality via ''m'' (mol/kg).', matA, matB);
    end
    rng_nacl = SF_phase_range('NaClaq');
    if any(m_val < rng_nacl.m(1)) || any(m_val > rng_nacl.m(2))
        error('SeaFreeze:badInput', ...
            'Molality value(s) outside the NaClaq spline range [%.3g, %.3g] mol/kg: %s.', ...
            rng_nacl.m(1), rng_nacl.m(2), mat2str(m_val));
    end
end

% ------------ Build the (P, T) sampling grid --------------------------------
rngA = SF_phase_range(matA);
rngB = SF_phase_range(matB);
P_lim = [max([rngA.P(1), rngB.P(1), 0.1]), min(rngA.P(2), rngB.P(2))];
T_lim = [max([rngA.T(1), rngB.T(1), 1.0]), min(rngA.T(2), rngB.T(2))];
if P_lim(1) >= P_lim(2) || T_lim(1) >= T_lim(2)
    error('SeaFreeze:unsupportedPair', ...
        'Phases %s and %s have non-overlapping spline domains.', matA, matB);
end

if isempty(P_user)
    dP = pick_step(P_lim(2) - P_lim(1), [500 1000 5000], [1 2 5 20]);
    P  = P_lim(1):dP:P_lim(2);
else
    P = P_user(:)';
end
if isempty(T_user)
    T = T_lim(1):0.5:T_lim(2);
else
    T = T_user(:)';
end

% ------------ Loop over molality values (single iter for pure pairs) --------
if nacl_involved
    m_iter = m_val;
else
    m_iter = NaN;     % single iteration with no molality
end
n_m = numel(m_iter);
results = [];   % initialised on first iteration so field set matches `r`

for k = 1:n_m
    [Ga, Gb] = compute_surfaces(matA, matB, P, T, m_iter(k));
    Z = Ga - Gb;
    TP = contourc(T, P, Z, [0 0]);
    [T_eq, P_eq] = parse_contourc(TP);

    if isempty(T_eq)
        if nacl_involved
            error('SeaFreeze:noContour', ...
                'No zero-crossing of G(%s)-G(%s) at m=%.3g mol/kg in the sampled grid.', ...
                matA, matB, m_iter(k));
        else
            error('SeaFreeze:noContour', ...
                'No zero-crossing of G(%s)-G(%s) found in the sampled (P, T) grid.', ...
                matA, matB);
        end
    end

    stable_mask = classify_stable(P_eq, T_eq, pair.stable_range);
    switch segment
        case 'all',    keep = true(size(stable_mask));
        case 'stable', keep = stable_mask;
        case 'meta',   keep = ~stable_mask;
    end

    r = struct();
    r.matA          = matA;
    r.matB          = matB;
    r.P             = P_eq(keep);
    r.T             = T_eq(keep);
    r.stable        = stable_mask(keep);
    r.segment       = segment;
    if nacl_involved, r.m = m_iter(k); end
    r.triple_points = pair.triple_points;
    if isempty(results)
        results = r;
    else
        results(k) = r;
    end
end

% ------------ Optional plot -------------------------------------------------
if do_plot
    fig = render_plot(results, target_fig);
    for k = 1:numel(results), results(k).fig = fig; end
end

% ------------ Return scalar struct for single curve, struct array for multi -
if numel(results) == 1
    out = results(1);
else
    out = results;
end
end


% ============================================================================
% Helpers
% ============================================================================

function step = pick_step(span, breaks, steps)
    % Adaptive step size: returns steps(k) for the first breaks(k) >= span.
    for k = 1:length(breaks)
        if span <= breaks(k), step = steps(k); return; end
    end
    step = steps(end);
end


function [Ga, Gb] = compute_surfaces(matA, matB, P, T, m)
    defs = sf_material_defs();
    MW_H2O = defs.MW_H2O;

    Ga = phase_surface(matA, P, T, m, MW_H2O);
    Gb = phase_surface(matB, P, T, m, MW_H2O);
end


function G = phase_surface(material, P, T, m, MW_H2O)
    % Returns a P-by-T matrix in J/mol-of-H2O.
    if strcmp(material, 'NaClaq')
        s = SF_getprop({P, T, m}, 'NaClaq', 'muw');
        % muw is P-by-T-by-1 because m is scalar; squeeze to 2D.
        G = squeeze(s.muw);
        if size(G, 1) == 1, G = G.'; end   % keep P down columns
    else
        s = SF_getprop({P, T}, material, 'G');
        G = s.G * MW_H2O;
    end
end


function [T_eq, P_eq] = parse_contourc(TP)
    % contourc output can chain multiple segments; each has a header
    % [level; npts]. Concatenate all (P, T) points across segments.
    T_eq = []; P_eq = [];
    if isempty(TP), return; end
    i = 1;
    while i <= size(TP, 2)
        npts = TP(2, i);
        T_eq = [T_eq; TP(1, i+1:i+npts).']; %#ok<AGROW>
        P_eq = [P_eq; TP(2, i+1:i+npts).']; %#ok<AGROW>
        i = i + npts + 1;
    end
end


function tf = valid_plot_arg(x)
    % 'plot' accepts: logical / scalar (true|false) or a figure handle.
    tf = islogical(x) || (isnumeric(x) && isscalar(x)) || ...
         isgraphics(x, 'figure');
end


function [Px, Tx] = group_with_breaks(P, T, mask)
    % Return the (P, T) values where mask==true, with NaN inserted between
    % each contiguous run. Lets a single plot() call render disjoint
    % segments (e.g. metastable extensions on either side of a stable
    % region) without drawing a fake connector across the gap.
    Px = []; Tx = [];
    if isempty(P) || ~any(mask), return; end
    m = logical(mask(:))';
    starts = find(diff([false m]) == 1);
    ends   = find(diff([m false]) == -1);
    P = P(:); T = T(:);
    for k = 1:length(starts)
        Px = [Px; P(starts(k):ends(k)); NaN]; %#ok<AGROW>
        Tx = [Tx; T(starts(k):ends(k)); NaN]; %#ok<AGROW>
    end
    if ~isempty(Px), Px(end) = []; Tx(end) = []; end   % strip trailing NaN
end


function stable = classify_stable(P_eq, T_eq, range)
    if isempty(range) || ~isstruct(range)
        stable = true(size(P_eq));
        return;
    end
    if range.lo > range.hi
        % Inverted bounds (e.g. lo=Inf, hi=-Inf) = "all metastable", used
        % for pairs whose contour is entirely metastable in the canonical
        % phase diagram (II ↔ water1). This branch must come before the
        % "all-non-finite" check below because Inf/-Inf are both non-finite.
        stable = false(size(P_eq));
        return;
    end
    if ~isfinite(range.lo) && ~isfinite(range.hi)
        % Both bounds non-finite (e.g. -Inf/Inf): treat whole curve as stable.
        % Used by NaClaq pairs where m-dependent triple points are out of scope.
        stable = true(size(P_eq));
        return;
    end
    if strcmp(range.var, 'T')
        stable = T_eq >= range.lo & T_eq <= range.hi;
    else
        stable = P_eq >= range.lo & P_eq <= range.hi;
    end
end


function fig = render_plot(results, target_fig)
    % Accepts either a scalar struct (single curve) or a struct array
    % (one entry per molality). target_fig is empty (= create new figure)
    % or an existing figure handle to overlay onto.
    if nargin < 2, target_fig = []; end

    if isempty(target_fig)
        fig = figure;
        ax  = gca(fig);
        new_axes = true;
    else
        fig = target_fig;
        ax  = gca(fig);
        new_axes = false;
    end
    hold(ax, 'on');
    n = numel(results);

    % Always set DisplayName on plotted lines so caller can build their own
    % legend across multiple SF_PhaseLines calls if they want.
    if n == 1
        out = results(1);
        if isfield(out,'m')
            label_base = sprintf('%s–%s (m=%.3g)', out.matA, out.matB, out.m);
        else
            label_base = sprintf('%s–%s', out.matA, out.matB);
        end
        if any(out.stable)
            [Ps, Ts] = group_with_breaks(out.P, out.T, out.stable);
            plot(ax, Ps, Ts, '-r', 'LineWidth', 1.4, 'DisplayName', label_base);
        end
        if any(~out.stable)
            [Pm, Tm] = group_with_breaks(out.P, out.T, ~out.stable);
            plot(ax, Pm, Tm, ':r', 'LineWidth', 1.0, ...
                 'DisplayName', [label_base ' (meta)'], 'HandleVisibility', 'off');
        end
        show_TPs = ~isfield(out,'m') || out.m <= 0.05;
        if show_TPs && ~isempty(out.triple_points)
            plot(ax, out.triple_points(:,2), out.triple_points(:,1), 'bo', ...
                 'MarkerFaceColor','b', 'MarkerSize', 5, 'HandleVisibility','off');
        end
        if new_axes
            title(ax, label_base);
        end
    else
        colors = lines(n);
        legend_handles = gobjects(0);
        legend_labels  = {};
        for k = 1:n
            r = results(k);
            label = sprintf('%s–%s (m=%.3g mol/kg)', r.matA, r.matB, r.m);
            h_stable = []; h_meta = [];
            if any(r.stable)
                [Ps, Ts] = group_with_breaks(r.P, r.T, r.stable);
                h_stable = plot(ax, Ps, Ts, '-', 'Color', colors(k,:), ...
                                'LineWidth', 1.4, 'DisplayName', label);
            end
            if any(~r.stable)
                [Pm, Tm] = group_with_breaks(r.P, r.T, ~r.stable);
                h_meta = plot(ax, Pm, Tm, ':', 'Color', colors(k,:), ...
                              'LineWidth', 1.0, 'HandleVisibility','off');
            end
            if ~isempty(h_stable)
                legend_handles(end+1) = h_stable; %#ok<AGROW>
            elseif ~isempty(h_meta)
                legend_handles(end+1) = h_meta; %#ok<AGROW>
            end
            legend_labels{end+1} = label; %#ok<AGROW>
        end
        all_low = all(arrayfun(@(r) ~isfield(r,'m') || r.m <= 0.05, results));
        if all_low && ~isempty(results(1).triple_points)
            plot(ax, results(1).triple_points(:,2), results(1).triple_points(:,1), ...
                 'bo', 'MarkerFaceColor','b', 'MarkerSize', 5, 'HandleVisibility','off');
        end
        if new_axes
            legend(ax, legend_handles, legend_labels, 'Location','best', 'FontSize', 9);
            title(ax, sprintf('%s — %s (%d molalities)', ...
                              results(1).matA, results(1).matB, n));
        end
    end

    if new_axes
        xlabel(ax, 'Pressure (MPa)');
        ylabel(ax, 'Temperature (K)');
        grid(ax, 'on');
    end
end


% ----------------------------------------------------------------------------
% Pair lookup table — single source of truth for triple points and
% stable-range bounds. Stable-range bounds are extracted from v1's magic-index
% substitutions (see SF_PhaseLines_v1.m). Pairs are symmetric.
% ----------------------------------------------------------------------------
function [pair, swapped] = lookup_pair(matA, matB)
    % Triple-point coordinates [T_K, P_MPa] (literature values used by v1).
    TP_IhLiqIII  = [251.165, 207.593];
    TP_IhIIIII   = [238.237, 209.885];
    TP_IIIIIV    = [249.418, 355.504];
    TP_IIVVI     = [201.934, 670.840];
    TP_IIIVLiq   = [256.164, 350.110];
    TP_VVILiq    = [273.407, 634.400];
    TP_VIVIILiq  = [354.000, 2216.000];   % VI–VII–liquid (Bridgman 1937 / IAPWS R14-08)
    TP_atm       = [273.150, 0.000611];   % atmospheric ice-Ih melting

    % Per-pair definitions: stable range as bounds in T (or P), and the
    % triple points to mark on plots.
    PAIRS = {
       %  matA      matB        var lo_value           hi_value           triple_points
        {'Ih',     'water1',    'T', TP_IhLiqIII(1),   TP_atm(1),         [TP_IhLiqIII; TP_atm]}
        {'Ih',     'II',        'T', 100.0,            TP_IhIIIII(1),     [TP_IhIIIII]}
        {'Ih',     'III',       'T', TP_IhIIIII(1),    TP_IhLiqIII(1),    [TP_IhIIIII; TP_IhLiqIII]}
        {'II',     'III',       'T', TP_IhIIIII(1),    TP_IIIIIV(1),      [TP_IhIIIII; TP_IIIIIV]}
        {'II',     'V',         'T', TP_IIVVI(1),      TP_IIIIIV(1),      [TP_IIVVI; TP_IIIIIV]}
        {'II',     'VI',        'T',  50.0,            TP_IIVVI(1),       [TP_IIVVI]}
        {'III',    'V',         'T', TP_IIIIIV(1),     TP_IIIVLiq(1),     [TP_IIIIIV; TP_IIIVLiq]}
        {'III',    'water1',    'T', TP_IhLiqIII(1),   TP_IIIVLiq(1),     [TP_IhLiqIII; TP_IIIVLiq]}
        {'V',      'water1',    'T', TP_IIIVLiq(1),    TP_VVILiq(1),      [TP_IIIVLiq; TP_VVILiq]}
        {'VI',     'water1',    'T', TP_VVILiq(1),     1000.0,            [TP_VVILiq]}
        {'V',      'VI',        'T', TP_IIVVI(1),      TP_VVILiq(1),      [TP_IIVVI; TP_VVILiq]}
        % II ↔ water1 — entirely metastable (II is surrounded by Ih/III/V/VI in
        % the standard phase diagram and never directly equilibrates with the
        % liquid). lo > hi triggers the "all-meta" branch in classify_stable.
        {'II',     'water1',    'T', Inf,              -Inf,              [TP_IhIIIII; TP_IIIIIV]}
        % VII_X_French ↔ liquid pairs. Stable above the VI–VII–water TP at
        % T = 354 K, P = 2216 MPa.
        {'VII_X_French','water1',         'P', TP_VIVIILiq(2), Inf,        [TP_VIVIILiq]}
        {'VII_X_French','water2',         'P', TP_VIVIILiq(2), Inf,        [TP_VIVIILiq]}
        {'VII_X_French','water_IAPWS95',  'P', TP_VIVIILiq(2), Inf,        [TP_VIVIILiq]}
        % Cross-EOS: same physical stable ranges as the water1 versions.
        % Useful for comparing alternative liquid parametrizations; results
        % drift from water1 at low P because water2 / IAPWS95 are not
        % optimised for the cold low-P regime.
        {'Ih',     'water2',         'T', TP_IhLiqIII(1),   TP_atm(1),     [TP_IhLiqIII; TP_atm]}
        {'Ih',     'water_IAPWS95',  'T', TP_IhLiqIII(1),   TP_atm(1),     [TP_IhLiqIII; TP_atm]}
        {'III',    'water2',         'T', TP_IhLiqIII(1),   TP_IIIVLiq(1), [TP_IhLiqIII; TP_IIIVLiq]}
        {'III',    'water_IAPWS95',  'T', TP_IhLiqIII(1),   TP_IIIVLiq(1), [TP_IhLiqIII; TP_IIIVLiq]}
        {'V',      'water2',         'T', TP_IIIVLiq(1),    TP_VVILiq(1),  [TP_IIIVLiq; TP_VVILiq]}
        {'V',      'water_IAPWS95',  'T', TP_IIIVLiq(1),    TP_VVILiq(1),  [TP_IIIVLiq; TP_VVILiq]}
        {'VI',     'water2',         'T', TP_VVILiq(1),     1000.0,        [TP_VVILiq]}
        {'VI',     'water_IAPWS95',  'T', TP_VVILiq(1),     1000.0,        [TP_VVILiq]}
        % NaClaq pairs — full curve marked stable by default. Distinguishing
        % stable vs metastable for ice ↔ NaClaq requires triple points that
        % depend on molality (out of scope for this rewrite). Triple-point
        % markers are inherited from the pure-ice ↔ water1 pair as a guide.
        {'Ih',     'NaClaq',    'T', -Inf,             Inf,               [TP_IhLiqIII; TP_atm]}
        {'III',    'NaClaq',    'T', -Inf,             Inf,               [TP_IhLiqIII; TP_IIIVLiq]}
        {'V',      'NaClaq',    'T', -Inf,             Inf,               [TP_IIIVLiq; TP_VVILiq]}
        {'VI',     'NaClaq',    'T', -Inf,             Inf,               [TP_VVILiq]}
        {'II',     'NaClaq',    'T', -Inf,             Inf,               [TP_IhIIIII]}
    };

    swapped = false;
    for k = 1:size(PAIRS, 1)
        row = PAIRS{k};
        if strcmp(matA, row{1}) && strcmp(matB, row{2})
            pair = make_pair_struct(row);
            return;
        elseif strcmp(matA, row{2}) && strcmp(matB, row{1})
            pair = make_pair_struct(row);
            swapped = true;
            return;
        end
    end
    error('SeaFreeze:unsupportedPair', ...
        'No equilibrium line defined for the pair (%s, %s).', matA, matB);
end


function pair = make_pair_struct(row)
    pair.matA          = row{1};
    pair.matB          = row{2};
    pair.stable_range  = struct('var', row{3}, 'lo', row{4}, 'hi', row{5});
    pair.triple_points = row{6};
end
