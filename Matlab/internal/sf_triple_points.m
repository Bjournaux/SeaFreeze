function tp = sf_triple_points(pm, fluid)
% SF_TRIPLE_POINTS  Triple points of a phase map (sf_phase_map), refined by
% Newton on G_a = G_b = G_c with dG/dP = 1/rho, dG/dT = -S.
%
% Candidates are the 2x2 grid blocks of pm.stable holding three distinct
% phases.  The ice Ih - liquid - vapour point (the fluid counts once in the
% map) is added where the saturation and sublimation curves cross.
%
% Output struct array: phases {a,b,c}, labels {..}, P (MPa), T (K),
%   rho [rho_a rho_b rho_c] (kg/m^3).
S = pm.stable;
B = cat(3, S(1:end-1, 1:end-1), S(2:end, 1:end-1), S(1:end-1, 2:end), S(2:end, 2:end));
srt = sort(B, 3);
nd = 1 + sum(diff(srt, 1, 3) ~= 0, 3);
[ii, jj] = find(nd >= 3 & srt(:, :, 1) > 0);
rows = zeros(numel(ii), 5);
for q = 1:numel(ii)
    u = unique(squeeze(srt(ii(q), jj(q), :)));
    rows(q, :) = [u(1:3).' log(pm.P(ii(q))) pm.T(jj(q))];
end
tp = struct('phases', {}, 'labels', {}, 'P', {}, 'T', {}, 'rho', {});
if ~isempty(rows)
    [keys, ~, g] = unique(rows(:, 1:3), 'rows');
    for k = 1:size(keys, 1)
        m = g == k;
        P = exp(median(rows(m, 4))); T = median(rows(m, 5));
        names = pm.names(keys(k, :));
        ok = false;
        for it = 1:30
            v = zeros(3, 3);
            for c = 1:3
                o = SF_getprop([P T], names{c}, {'G', 'S', 'rho'});
                v(c, :) = [o.G o.S o.rho];
                if o.G == 0, v(c, :) = NaN; end          % outside a Gibbs spline
            end
            if any(~isfinite(v(:))), break; end
            r = [v(1,1) - v(2,1); v(1,1) - v(3,1)];
            J = [1e6 * (1/v(1,3) - 1/v(2,3)), -(v(1,2) - v(2,2));
                 1e6 * (1/v(1,3) - 1/v(3,3)), -(v(1,2) - v(3,2))];
            d = -J \ r;
            if any(~isfinite(d)), break; end
            dP = max(min(d(1), 0.5 * P), -0.5 * P); dT = max(min(d(2), 10), -10);
            P = P + dP; T = T + dT;
            if abs(dP) < 1e-9 * max(P, 1e-6) && abs(dT) < 1e-9, ok = true; break; end
        end
        if ok
            rho = zeros(1, 3);
            for c = 1:3
                o = SF_getprop([P T], names{c}, 'rho'); rho(c) = o.rho;
            end
            lab = strrep(names, 'VII_X_French', 'VII/X');
            if keys(k, 1) == 1, lab{1} = 'L'; end
            tp(end+1) = struct('phases', {names}, 'labels', {lab}, 'P', P, 'T', T, 'rho', rho); %#ok<AGROW>
        end
    end
end
% ice Ih - liquid - vapour: saturation meets sublimation
if any(strcmp(pm.names, 'Ih')) && min(pm.T) < 273.16 && max(pm.T) > 273.16
    f = @(T) log(getfield(SF_coexistence('saturation', T, 'fluid', fluid), 'P')) - ...
             log(getfield(SF_coexistence('sublimation', T, 'fluid', fluid), 'P'));
    a = 265; b = 280; fa = f(a);
    for it = 1:50
        mdl = 0.5 * (a + b); fm = f(mdl);
        if sign(fm) == sign(fa), a = mdl; fa = fm; else, b = mdl; end
    end
    Tt = 0.5 * (a + b);
    st = SF_coexistence('saturation', Tt, 'fluid', fluid);
    oI = SF_getprop([st.P Tt], 'Ih', 'rho');
    tp(end+1) = struct('phases', {{'Ih', fluid, fluid}}, 'labels', {{'Ih', 'L', 'V'}}, ...
                       'P', st.P, 'T', Tt, 'rho', [oI.rho st.rho_A st.rho_B]);
end
end
