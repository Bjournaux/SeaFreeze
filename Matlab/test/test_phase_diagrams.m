function test_phase_diagrams()
% Tests for sf_phase_map, sf_triple_points, SF_WPD_PT and SF_WPD_rhoT.
% Baptiste Journaux - 2026
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'internal'));
addpath(fullfile(fileparts(here), 'internal', 'compat'));
addpath(fileparts(here));
warning('off', 'SeaFreeze:longRuntime');
np = 0; nf = 0;

% spot checks of the stable phase
m = sf_phase_map([1e-6 0.1 0.1 300 1000], [300 260 300 230 250]);
st = arrayfun(@(i) m.stable(i, i), 1:5);
[np,nf] = check('phase map spot checks (vapour, Ih, liquid, II, VI)', ...
    isequal(m.names(st), {'water_Brown2026','Ih','water_Brown2026','II','VI'}) && m.rho_stable(1,1) < 1e-3 && m.rho_stable(3,3) > 990, np, nf);

% triple points vs literature (SF_PhaseLines table)
pm = sf_phase_map(logspace(-9, log10(3000), 500), linspace(180, 420, 240));
tp = sf_triple_points(pm, 'water_Brown2026');
lit = {{'Ih','II','III'}, 238.237, 209.885; {'II','III','V'}, 249.418, 355.504; ...
       {'II','V','VI'}, 201.934, 670.840; {'Ih','III','L'}, 251.165, 207.593; ...
       {'III','L','V'}, 256.164, 350.110; {'L','V','VI'}, 273.407, 634.400};
ok = true;
for k = 1:size(lit, 1)
    hit = false;
    for q = 1:numel(tp)
        if isequal(sort(tp(q).labels), sort(lit{k,1}))
            hit = abs(tp(q).T - lit{k,2}) < 0.1 && abs(tp(q).P - lit{k,3}) < 1.5;
            if hit, break; end
        end
    end
    ok = ok && hit;
end
lv = tp(arrayfun(@(t) isequal(sort(t.labels), {'Ih','L','V'}), tp));
[np,nf] = check(sprintf('triple points vs literature (%d found)', numel(tp)), ok, np, nf);
[np,nf] = check(sprintf('Ih-L-V triple point %.4f K, %.3f Pa', lv.T, lv.P * 1e6), ...
    numel(lv) == 1 && abs(lv.T - 273.16) < 0.01 && abs(lv.P * 1e6 - 611.657) < 0.5, np, nf);

% the diagram functions run (small grids)
f1 = SF_WPD_PT('nP', 120, 'nT', 90); set(f1, 'Visible', 'off');
f2 = SF_WPD_rhoT('nP', 240, 'nT', 90, 'nrho', 200); set(f2, 'Visible', 'off');
f3 = SF_WPD_rhoT('rho', [850 1700], 'T', [150 500], 'P', [1e-10 1e4], 'xscale', 'linear', ...
                 'nP', 300, 'nT', 90, 'nrho', 200); set(f3, 'Visible', 'off');
[np,nf] = check('SF_WPD_PT / SF_WPD_rhoT render', all(ishghandle([f1 f2 f3])), np, nf);
close([f1 f2 f3]);

% melting curve (psiEOS dq2026)
Tm = sf_melt_T_dq2026([1e-4 0.101325 300 1000 10000]);
[np,nf] = check('sf_melt_T_dq2026', abs(Tm(1) - 273.16) < 1e-9 && abs(Tm(2) - 273.152) < 0.01 && ...
    Tm(3) > 250 && Tm(3) < 260 && Tm(4) > 300 && Tm(4) < 310 && Tm(5) > 600 && Tm(5) < 800, np, nf);

warning('on', 'SeaFreeze:longRuntime');
fprintf('\n%d passed, %d failed\n', np, nf);
if nf > 0, error('test_phase_diagrams:failed', '%d test(s) failed', nf); end
end

function [np, nf] = check(name, cond, np, nf)
if cond, fprintf('[PASS] %s\n', name); np = np + 1;
else,    fprintf('[FAIL] %s\n', name); nf = nf + 1; end
end
