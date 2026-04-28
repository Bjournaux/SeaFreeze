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
addParameter(p, 'm',       NaN,    @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'P',       [],     @(x) isnumeric(x) && isvector(x));
addParameter(p, 'T',       [],     @(x) isnumeric(x) && isvector(x));
addParameter(p, 'segment', 'all',  @(s) ismember(lower(char(s)), {'all','stable','meta'}));
addParameter(p, 'plot',    false,  @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
parse(p, matA, matB, varargin{:});
matA = char(matA);
matB = char(matB);
m_val   = p.Results.m;
P_user  = p.Results.P;
T_user  = p.Results.T;
segment = lower(char(p.Results.segment));
do_plot = logical(p.Results.plot);

% ------------ Look up the pair (handles symmetric ordering) -----------------
[pair, swapped] = lookup_pair(matA, matB);
if swapped, [matA, matB] = deal(matB, matA); end

% ------------ Validate molality if NaClaq is involved -----------------------
nacl_involved = any(strcmp({matA, matB}, 'NaClaq'));
if nacl_involved
    if isnan(m_val)
        error('SeaFreeze:badInput', ...
            'Pair (%s, %s) involves NaClaq; pass molality via ''m'' (mol/kg).', matA, matB);
    end
    rng_nacl = SF_phase_range('NaClaq');
    if m_val < rng_nacl.m(1) || m_val > rng_nacl.m(2)
        error('SeaFreeze:badInput', ...
            'Molality %.3g is outside the NaClaq spline range [%.3g, %.3g] mol/kg.', ...
            m_val, rng_nacl.m(1), rng_nacl.m(2));
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

% ------------ Compute G surfaces and zero contour ---------------------------
[Ga, Gb] = compute_surfaces(matA, matB, P, T, m_val);
Z = Ga - Gb;
TP = contourc(T, P, Z, [0 0]);
% contourc may return multiple disjoint segments concatenated, each with a
% header column [level; npts]. Flatten and skip the headers.
[T_eq, P_eq] = parse_contourc(TP);

if isempty(T_eq)
    error('SeaFreeze:noContour', ...
        'No zero-crossing of G(%s) - G(%s) found in the sampled (P, T) grid.', matA, matB);
end

% ------------ Stable / metastable classification ----------------------------
stable_mask = classify_stable(P_eq, T_eq, pair.stable_range);

% ------------ Apply segment filter ------------------------------------------
switch segment
    case 'all',    keep = true(size(stable_mask));
    case 'stable', keep = stable_mask;
    case 'meta',   keep = ~stable_mask;
end

% ------------ Build output struct -------------------------------------------
out = struct();
out.matA    = matA;
out.matB    = matB;
out.P       = P_eq(keep);
out.T       = T_eq(keep);
out.stable  = stable_mask(keep);
out.segment = segment;
if nacl_involved, out.m = m_val; end
out.triple_points = pair.triple_points;   % Nx2 array, columns [T_K, P_MPa]

% ------------ Optional plot -------------------------------------------------
if do_plot
    out.fig = render_plot(out);
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
    MW_H2O = 0.018015268;     % kg/mol; matches fnGval / SF_WhichPhase

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


function stable = classify_stable(P_eq, T_eq, range)
    if isempty(range) || ~isstruct(range) || ...
            (~isfinite(range.lo) && ~isfinite(range.hi))
        % No bounds defined: treat the whole curve as stable.
        stable = true(size(P_eq));
        return;
    end
    if strcmp(range.var, 'T')
        stable = T_eq >= range.lo & T_eq <= range.hi;
    else
        stable = P_eq >= range.lo & P_eq <= range.hi;
    end
end


function fig = render_plot(out)
    fig = figure;
    hold on;
    if any(out.stable)
        plot(out.P(out.stable), out.T(out.stable), '-r', 'LineWidth', 1.4);
    end
    if any(~out.stable)
        plot(out.P(~out.stable), out.T(~out.stable), ':r', 'LineWidth', 1.0);
    end
    if ~isempty(out.triple_points)
        plot(out.triple_points(:,2), out.triple_points(:,1), 'bo', ...
             'MarkerFaceColor', 'b', 'MarkerSize', 5);
    end
    xlabel('Pressure (MPa)');
    ylabel('Temperature (K)');
    if isfield(out, 'm')
        title(sprintf('%s — %s (m = %.3g mol/kg)', out.matA, out.matB, out.m));
    else
        title(sprintf('%s — %s', out.matA, out.matB));
    end
    grid on;
    hold off;
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
