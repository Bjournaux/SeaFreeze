function sf_wpd_labels(ax, X, Y, pm, kind, img)
% SF_WPD_LABELS  Field labels at the median position of each stability field.
%   kind 'PT'  : X = log10(P), Y = T, fields from pm.stable
%   kind 'rhoT': X = plotted density coordinate, Y = T, fields from img
Tc = 647.096; rhoc = 322;
[XX, YY] = ndgrid(X(:), Y(:));
if strcmp(kind, 'PT')
    S = pm.stable;
    R = pm.rho_stable;
else
    S = img; R = [];
end
lab = strrep(pm.names, 'VII_X_French', 'VII/X');
for k = 2:numel(pm.names)
    m = S == k;
    if nnz(m) > 30
        text(ax, median(XX(m)), median(YY(m)), lab{k}, 'FontWeight', 'bold', 'FontSize', 10, ...
             'HorizontalAlignment', 'center');
    end
end
m0 = S == 1;
if strcmp(kind, 'PT')
    sel = {m0 & YY < Tc & R < rhoc, m0 & YY < Tc & R >= rhoc, m0 & YY >= Tc};
else
    sel = {m0 & YY < Tc & XX < pm.xrhoc, m0 & YY < Tc & XX >= pm.xrhoc, m0 & YY >= Tc + 50};
end
names = {'vapour', 'liquid', 'supercritical fluid'};
for k = 1:3
    if nnz(sel{k}) > 150
        text(ax, median(XX(sel{k})), median(YY(sel{k})), names{k}, 'FontWeight', 'bold', ...
             'FontSize', 11, 'Color', [0.07 0.24 0.44], 'HorizontalAlignment', 'center');
    end
end
if strcmp(kind, 'PT')
    mm = S == 0;
    if nnz(mm) > 30
        text(ax, median(XX(mm)), median(YY(mm)), {'not modelled', '(ice VII/X field)'}, ...
             'FontSize', 9, 'Color', [0.32 0.32 0.31], 'HorizontalAlignment', 'center');
    end
end
end
