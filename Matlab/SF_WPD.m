function fig = SF_WPD(varargin)
% SF_WPD  Draw the H2O Water Phase Diagram.
% Version 2.0.0
% Baptiste Journaux - 2026
%
% Computes all phase boundaries dynamically from Gibbs energy splines
% via SF_PhaseLines.  Supports NaCl(aq) melting-curve overlays.
%
% Usage:
%   SF_WPD()                                        % pure water, new figure
%   SF_WPD('solute','NaCl', 'm', [0.5 1 2 4])      % NaClaq overlay
%   SF_WPD('meta', false)                           % hide metastable extensions
%   SF_WPD('labels', true)                          % annotate phase fields
%   SF_WPD('ax', gca)                               % overlay on existing axes
%   fig = SF_WPD(...)                               % return figure handle
%
% Name-Value Parameters:
%   'ax'      - axes handle to plot onto (default: new figure)
%   'solute'  - 'none' (default) or 'NaCl' to overlay NaClaq melting curves
%   'm'       - scalar or vector of molality values (mol/kg) for NaClaq
%   'meta'    - 'default' (only Ih-II and II-VI, matching v1), true/'all'
%               (all pairs), or false/'none' (no metastable extensions)
%   'labels'  - logical, annotate stability fields with phase names (default: false)
%
% See also: SF_PhaseLines, SF_WhichPhase

p = inputParser;
addParameter(p, 'ax',      [],      @(x) isempty(x) || isgraphics(x, 'axes'));
addParameter(p, 'solute',  'none',  @(s) ischar(s) || (isstring(s) && isscalar(s)));
addParameter(p, 'm',       [],      @(x) isnumeric(x) && (isempty(x) || isvector(x)));
addParameter(p, 'meta',    'default', @valid_meta);
addParameter(p, 'labels',  true,    @(x) islogical(x) || (isnumeric(x) && isscalar(x)));
parse(p, varargin{:});

ax_in     = p.Results.ax;
solute    = lower(char(p.Results.solute));
m_arr     = p.Results.m(:)';
meta_arg  = p.Results.meta;
do_labels = logical(p.Results.labels);
is_nacl   = ismember(solute, {'nacl','naclaq'});

% Parse meta option: 'default' | true/'all' | false/'none'
DEFAULT_META_PAIRS = {'Ih','II'; 'II','VI'};
if ischar(meta_arg) || isstring(meta_arg)
    switch lower(char(meta_arg))
        case 'default', meta_mode = 1;   % only Ih-II and II-VI
        case 'all',     meta_mode = 2;   % all pairs
        case 'none',    meta_mode = 0;   % no metastable
        otherwise,      meta_mode = 1;
    end
elseif meta_arg
    meta_mode = 2;  % true  -> all
else
    meta_mode = 0;  % false -> none
end

if is_nacl && isempty(m_arr)
    error('SF_WPD:badInput', ...
          'solute is set to NaCl; pass molality value(s) via ''m''.');
end

% --- Create or reuse axes -----------------------------------------------------
if isempty(ax_in)
    fig = figure('Position', [100 100 900 650]);
    ax  = axes(fig);
else
    ax  = ax_in;
    fig = ancestor(ax, 'figure');
end
hold(ax, 'on');

% --- Pure-water / ice phase boundaries ----------------------------------------
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
    plot_runs(ax, r.P, r.T, r.stable, '-', 1.4, [0 0 0]);
    % Metastable: mode 2 = all pairs; mode 1 = only Ih-II and II-VI below 150 K
    if meta_mode == 2
        plot_runs(ax, r.P, r.T, ~r.stable, '--', 1.0, [0.5 0.5 0.5]);
    elseif meta_mode == 1 && is_default_meta(PAIRS{k,1}, PAIRS{k,2}, DEFAULT_META_PAIRS)
        low_T = ~r.stable & r.T < 150;
        plot_runs(ax, r.P, r.T, low_T, '--', 1.0, [0.5 0.5 0.5]);
    end
end

% --- NaClaq melting-curve overlay ---------------------------------------------
if is_nacl
    nm = numel(m_arr);
    cmap = parula(max(nm, 2));
    ice_phases = {'Ih', 'III', 'V', 'VI'};
    for j = 1:nm
        for ic = 1:numel(ice_phases)
            try
                r = SF_PhaseLines(ice_phases{ic}, 'NaClaq', 'm', m_arr(j));
            catch
                continue
            end
            lbl = '';
            if ic == 1
                lbl = sprintf('m = %g mol/kg', m_arr(j));
            end
            plot(ax, r.P, r.T, '-', 'Color', cmap(j,:), ...
                 'LineWidth', 1.2, 'DisplayName', lbl);
        end
    end
    legend(ax, 'show', 'Location', 'northeast', 'FontSize', 9);
end

% --- Axis formatting ----------------------------------------------------------
xlabel(ax, 'Pressure (MPa)', 'FontSize', 12);
ylabel(ax, 'Temperature (K)', 'FontSize', 12);
title(ax, 'H_2O Phase Diagram (SeaFreeze)', 'FontSize', 13);
xlim(ax, [0 2300]);
ylim(ax, [100 400]);
grid(ax, 'on');
set(ax, 'GridAlpha', 0.3);

% --- Phase-field labels -------------------------------------------------------
if do_labels
    kw = {'FontSize', 11, 'FontWeight', 'bold', ...
          'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
          'Color', [0.2 0.2 0.2], ...
          'BackgroundColor', 'w', 'EdgeColor', 'none', 'Margin', 2};
    text(ax, 80,   240, 'Ih',     kw{:});
    text(ax, 330,  217, 'II',     kw{:});
    text(ax, 280,  248, 'III',    kw{:});
    text(ax, 490,  247, 'V',      kw{:});
    text(ax, 1200, 255, 'VI',     kw{:});
    text(ax, 300,  345, 'Liquid', kw{:});
end

if nargout == 0
    clear fig
end

end


% ==============================================================================
function tf = valid_meta(x)
    if islogical(x) || (isnumeric(x) && isscalar(x))
        tf = true;
    elseif ischar(x) || (isstring(x) && isscalar(x))
        tf = ismember(lower(char(x)), {'default','all','none'});
    else
        tf = false;
    end
end

function tf = is_default_meta(a, b, pairs)
    tf = false;
    for j = 1:size(pairs, 1)
        if (strcmp(a, pairs{j,1}) && strcmp(b, pairs{j,2})) || ...
           (strcmp(a, pairs{j,2}) && strcmp(b, pairs{j,1}))
            tf = true; return;
        end
    end
end

function plot_runs(ax, P, T, mask, style, lw, color)
% Insert NaN between disjoint runs so a single plot() call doesn't bridge gaps.
    idx = find(mask);
    if isempty(idx), return; end
    breaks = find(diff(idx) > 1);
    Pp = P(idx);
    Tp = T(idx);
    for b = flip(breaks(:)')
        Pp = [Pp(1:b); NaN; Pp(b+1:end)];
        Tp = [Tp(1:b); NaN; Tp(b+1:end)];
    end
    plot(ax, Pp, Tp, style, 'LineWidth', lw, 'Color', color, ...
         'HandleVisibility', 'off');
end
