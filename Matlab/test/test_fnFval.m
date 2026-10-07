function test_fnFval()
% Tests for fnFval / psi_val — Helmholtz-spline water (water_Brown2026) and the plain
% F(rho,T) test spline.
% Baptiste Journaux - 2026
%
% Fixtures (test/fixtures/):
%   psi_reference.mat   psiH2O_val.m (lbf-thermo, its own engine/lbf/sp_val) output
%                       at scattered / gridded / edge states of the water_Brown2026 surface
%   water_F_test.mat    F(rho,T) B-spline fitted to IAPWS-95 (gen_helmholtz_fixture.m)

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'internal'));
addpath(fullfile(fileparts(here), 'internal', 'compat'));
addpath(fileparts(here));

np = 0; nf = 0;
relerr = @(a, b) max(abs(a(:) - b(:)) ./ max(abs(b(:)), 1e-300));

% =========================================================================
% 1. psi_val vs psiH2O_val reference (rho,T mode)
% =========================================================================
R = load(fullfile(here, 'fixtures', 'psi_reference.mat')); ref = R.ref;
sp = sf_load_spline('water_Brown2026');
[np,nf] = check('water_Brown2026 spline carries eos = psi', isfield(sp,'eos') && strcmp(sp.eos,'psi'), np, nf);

tic;
o = fnFval(sp, [ref.rho ref.T], [], 'rhoT');
t1 = toc;
fprintf('   (fnFval rhoT scatter, %d states, all props: %.2f s)\n', numel(ref.rho), t1);
r = ref.out;
tol = 1e-9;
pairs = {'P','P'; 'G','G'; 'S','S'; 'U','U'; 'H','H'; 'A','F'; 'Cp','cp'; 'Cv','cv'; ...
         'Kt','Kt'; 'Ks','Ks'; 'alpha','alpha'};
for k = 1:size(pairs,1)
    a = o.(pairs{k,1}); b = real(r.(pairs{k,2}));
    ok = isfinite(b);
    e = relerr(a(ok), b(ok));
    [np,nf] = check(sprintf('psi scatter: %-5s rel err %.1e', pairs{k,1}, e), e < tol && isequal(isnan(a), isnan(b)), np, nf);
end
% sound speed: reference is complex where w^2 < 0
wv = r.w; okw = isfinite(wv) & imag(wv) == 0;
e = relerr(o.vel(okw), real(wv(okw)));
[np,nf] = check(sprintf('psi scatter: vel   rel err %.1e', e), e < tol, np, nf);
[np,nf] = check('psi scatter: vel NaN where w^2 < 0', all(isnan(o.vel(imag(wv) ~= 0))), np, nf);

% gridded rhoT
g = fnFval(sp, {ref.grid.rho, ref.grid.T}, {'P','G','Cp','rho','T'}, 'rhoT');
e = max(relerr(g.P, ref.grid.out.P), relerr(g.G, ref.grid.out.G));
[np,nf] = check(sprintf('psi grid: shape and values (rel err %.1e)', e), ...
    isequal(size(g.P), [numel(ref.grid.rho) numel(ref.grid.T)]) && e < tol && ...
    isequal(g.rho, ref.grid.rho(:)) && isequal(g.T, ref.grid.T(:)), np, nf);

% edge states: below floor (virial continuation) and outside (NaN)
q = fnFval(sp, [ref.edge.rho ref.edge.T], {'P','G'}, 'rhoT');
[np,nf] = check('psi edge: below-floor state matches reference', ...
    relerr(q.P(1), ref.edge.out.P(1)) < tol && relerr(q.G(1), ref.edge.out.G(1)) < tol, np, nf);
[np,nf] = check('psi edge: outside box -> NaN', all(isnan(q.P(2:4))) && all(isnan(q.G(2:4))), np, nf);

% Kp by finite difference: compare with FD of Kt in P
pt = [1000 300; 1200 350; 950 450];
o1 = fnFval(sp, pt, {'Kp','Kt','P'}, 'rhoT');
h = 1e-3;
op = fnFval(sp, [pt(:,1)*(1+h) pt(:,2)], {'Kt','P'}, 'rhoT');
om = fnFval(sp, [pt(:,1)*(1-h) pt(:,2)], {'Kt','P'}, 'rhoT');
Kp_fd = (op.Kt - om.Kt) ./ (op.P - om.P);
[np,nf] = check(sprintf('psi Kp consistent with dKt/dP (rel err %.1e)', relerr(o1.Kp, Kp_fd)), relerr(o1.Kp, Kp_fd) < 1e-4, np, nf);

% =========================================================================
% 2. PT mode: inversion round trip and SF_getprop dispatch
% =========================================================================
P = [0.1 1 10 100 500 1000 2000]; T = [240 273.16 300 350 400 500];
tic;
w = SF_getprop({P, T}, 'water_Brown2026');
t2 = toc;
fprintf('   (SF_getprop water_Brown2026 PT grid %dx%d, all props: %.2f s)\n', numel(P), numel(T), t2);
[np,nf] = check('SF_getprop water_Brown2026: grid shape', isequal(size(w.rho), [numel(P) numel(T)]) && ...
    isequal(w.P, P(:)) && isequal(w.T, T(:)), np, nf);
[Pm, Tm] = ndgrid(P, T);
fin = isfinite(w.rho(:));
fprintf('   (%d of %d grid points have a liquid-branch root)\n', nnz(fin), numel(fin));
if ~all(fin), disp([Pm(~fin) Tm(~fin)]); end
b = SF_getprop([w.rho(fin) Tm(fin)], 'water_Brown2026', 'P', 'input', 'rhoT');
[np,nf] = check(sprintf('PT -> rho -> P round trip (max |dP| %.1e MPa)', max(abs(b.P - Pm(fin)))), ...
    max(abs(b.P - Pm(fin))) < 1e-7, np, nf);
% liquid at ambient
a = SF_getprop([0.101325 298.15], 'water_Brown2026');
[np,nf] = check(sprintf('ambient liquid: rho = %.3f, Cp = %.1f', a.rho, a.Cp), ...
    abs(a.rho - 997.05) < 0.1 && abs(a.Cp - 4181.5) < 5, np, nf);
% thermodynamic identities
s = SF_getprop([50 300; 500 350; 1500 400; 0.5 260], 'water_Brown2026');
id1 = relerr(s.Cp - s.Cv, s.T .* s.alpha.^2 .* s.Kt * 1e6 ./ s.rho);
id2 = relerr(s.Ks ./ s.Kt, s.Cp ./ s.Cv);
id3 = relerr(s.vel.^2, s.Ks * 1e6 ./ s.rho);
id4 = relerr(s.G, s.U - s.T .* s.S + s.P * 1e6 ./ s.rho);
[np,nf] = check('identities Cp-Cv, Ks/Kt, vel, G', max([id1 id2 id3 id4]) < 1e-9, np, nf);
% scatter vs grid consistency
sc = SF_getprop([P(3) T(3)], 'water_Brown2026');
[np,nf] = check('scatter == grid at a shared point', abs(sc.rho - w.rho(3,3)) < 1e-9 && abs(sc.G - w.G(3,3)) < 1e-6, np, nf);
% out of range -> NaN, no error
z = SF_getprop([100 200; 1e7 300], 'water_Brown2026', {'rho','G'});
[np,nf] = check('out-of-domain PT -> NaN', all(isnan(z.rho)) && all(isnan(z.G)), np, nf);
% branch option refused for Gibbs materials
try
    SF_getprop([100 300], 'water_Bollengier2019', 'rho', 'branch', 'liquid'); bad = false;
catch e
    bad = strcmp(e.identifier, 'SeaFreeze:badInput');
end
[np,nf] = check('branch refused for Gibbs materials', bad, np, nf);

% (rho,T) input for Gibbs splines: P from SF_rho2P, then (P,T) evaluation
g1 = SF_getprop({[1000 1050 1100], [280 300 330]}, 'water_Bollengier2019', {'P','rho','G','Cp','T'}, 'input', 'rhoT');
[Pg, Tg] = deal(g1.P(:), reshape(repmat([280 300 330], 3, 1), [], 1));
chk1 = SF_getprop([Pg Tg], 'water_Bollengier2019', {'rho','G','Cp'});
[np,nf] = check(sprintf('rhoT grid water_Bollengier2019: shape, rho round trip %.1e', max(abs(chk1.rho - reshape(repmat([1000;1050;1100],1,3),[],1)))), ...
    isequal(size(g1.P), [3 3]) && isequal(g1.rho, [1000;1050;1100]) && isequal(g1.T, [280;300;330]) && ...
    max(abs(chk1.rho - reshape(repmat([1000;1050;1100],1,3),[],1))) < 1e-4 && ...
    max(abs(chk1.G - g1.G(:))) < 1e-9 * max(abs(g1.G(:))) + 1e-9, np, nf);
g6 = SF_getprop([1330 260; 1360 270; 2000 270], 'VI', {'P','Vp','shear'}, 'input', 'rhoT');
[np,nf] = check('rhoT scatter ice VI incl. shear/Vp; out-of-range -> NaN', ...
    all(isfinite(g6.P(1:2))) && all(isfinite(g6.Vp(1:2))) && isnan(g6.P(3)) && isnan(g6.Vp(3)), np, nf);
gn = SF_getprop([1050 300 1; 1100 320 2], 'NaClaq_Brown2026', {'P','rho','muw'}, 'input', 'rhoT');
bk = SF_getprop([gn.P gn.rho(:)*0 + [300; 320] [1; 2]], 'NaClaq_Brown2026', 'rho');
[np,nf] = check('rhoT scatter NaClaq (P,T,m) round trip', max(abs(bk.rho - [1050; 1100])) < 1e-4, np, nf);

% branch selection: vapour below the saturation pressure, liquid above
v  = SF_getprop([1e-3 300; 0.1 400; 0.1 300; 10 400], 'water_Brown2026', {'rho','G'});
vl = SF_getprop([1e-3 300; 0.1 400], 'water_Brown2026', {'rho','G'}, 'branch', 'liquid');
vv = SF_getprop([0.1 300], 'water_Brown2026', {'rho','G'}, 'branch', 'vapor');
[np,nf] = check(sprintf('stable branch: vapour at 1 kPa/300 K (rho %.4g) and 0.1 MPa/400 K (rho %.4g)', v.rho(1), v.rho(2)), ...
    v.rho(1) < 0.01 && v.rho(2) < 1 && v.rho(3) > 990 && v.rho(4) > 930, np, nf);
[np,nf] = check('liquid branch: metastable liquid with higher G', all(vl.rho > 930) && all(vl.G > v.G(1:2)), np, nf);
[np,nf] = check('vapor branch: metastable vapour with higher G at 0.1 MPa/300 K', vv.rho < 1 && vv.G > v.G(3), np, nf);

% =========================================================================
% 3. SF_rho2P / SF_phase_range for the Helmholtz material
% =========================================================================
Pr = SF_rho2P([997.05 1100 1200], [298.15 300 350], 'water_Brown2026');
chk = SF_getprop([Pr(:) [298.15;300;350]], 'water_Brown2026', 'rho');
[np,nf] = check('SF_rho2P water_Brown2026 round trip', max(abs(chk.rho(:) - [997.05;1100;1200])) < 1e-6, np, nf);
rg = SF_phase_range('water_Brown2026');
[np,nf] = check('SF_phase_range water_Brown2026 has rho, T, P', all(isfield(rg, {'rho','T','P'})) && rg.T(1) < 240 && rg.P(2) > 2300, np, nf);

% =========================================================================
% 4. Plain F(rho,T) spline (IAPWS-95 fit) — generic path
% =========================================================================
F = load(fullfile(here, 'fixtures', 'water_F_test.mat')); spF = F.sp;
o = fnFval(spF, {[1000 1100 1200], [280 300 400]}, [], 'rhoT');
[np,nf] = check('F_rhoT: finite output on grid', all(isfinite(o.P(:))) && all(isfinite(o.Cp(:))), np, nf);
q = fnFval(spF, {[1 100 1000], [280 300 400]}, {'rho','P'});
b = fnFval(spF, [q.rho(:) reshape(repmat([280 300 400], 3, 1), [], 1)], 'P', 'rhoT');
[np,nf] = check('F_rhoT: PT round trip', max(abs(b.P - reshape(repmat([1;100;1000], 1, 3), [], 1))) < 1e-7, np, nf);
% cross-check the two representations at a liquid state (both IAPWS-95-like)
a1 = fnFval(spF, [100 300], {'rho','Cp','S'}); a2 = SF_getprop([100 300], 'water_Brown2026', {'rho','Cp','S'});
[np,nf] = check(sprintf('F_rhoT fit vs water_Brown2026 at 100 MPa/300 K: drho %.2e', abs(a1.rho/a2.rho-1)), ...
    abs(a1.rho/a2.rho - 1) < 1e-3 && abs(a1.Cp/a2.Cp - 1) < 1e-2, np, nf);

% =========================================================================
% 5. Phase equilibria with water_Brown2026 as the liquid
% =========================================================================
PT = {[0.1 100 300 800 1500], [250 260 270 276 300 330]};
w3 = SF_WhichPhase(PT, 'liquid', 'water_Brown2026');
w1 = SF_WhichPhase(PT);
[np,nf] = check('SF_WhichPhase liquid water_Brown2026 == water_Bollengier2019 (grid)', isequal(w3, w1), np, nf);
ws = SF_WhichPhase([0.1 270; 0.1 276; 1000 300], 'liquid', 'water_Brown2026');
[np,nf] = check('SF_WhichPhase liquid water_Brown2026 (scatter)', isequal(ws(:).', [1 0 6]), np, nf);
r = SF_PhaseLines('Ih', 'water_Brown2026', 'P', [0.1 0.101325 1], 'T', 270:0.005:276);
[Ps, io] = sort(r.P); T01 = interp1(Ps, r.T(io), 0.101325);
[np,nf] = check(sprintf('Ih-water_Brown2026 melting at 0.101325 MPa: %.4f K', T01), abs(T01 - 273.152) < 0.01, np, nf);
for ice = {'Ih','III','V','VI'}
    a3 = SF_PhaseLines(ice{1}, 'water_Brown2026', 'segment', 'stable');
    a1 = SF_PhaseLines(ice{1}, 'water_Bollengier2019', 'segment', 'stable');
    [P3, i3] = sort(a3.P); [P1, i1] = sort(a1.P);
    [P3, u3] = unique(P3); T3 = a3.T(i3); T3 = T3(u3);
    [P1, u1] = unique(P1); T1 = a1.T(i1); T1 = T1(u1);
    Pc = linspace(max(P3(1), P1(1)), min(P3(end), P1(end)), 30);
    dT = max(abs(interp1(P3, T3, Pc) - interp1(P1, T1, Pc)));
    tol = 0.03; if strcmp(ice{1}, 'V'), tol = 0.08; elseif strcmp(ice{1}, 'VI'), tol = 0.5; end
    [np,nf] = check(sprintf('%s-liquid line water_Brown2026 vs water_Bollengier2019: max |dT| %.3f K', ice{1}, dT), dT < tol, np, nf);
end

% =========================================================================
% 6. Saturation and sublimation (SF_coexistence)
% =========================================================================
Ts = [273.16 300 373.124 450 550 620];
sa = SF_coexistence('saturation', Ts);
aux = 22.064 * exp(647.096 ./ Ts(:) .* (-7.85951783*(1-Ts(:)/647.096) + 1.84408259*(1-Ts(:)/647.096).^1.5 ...
      - 11.7866497*(1-Ts(:)/647.096).^3 + 22.6807411*(1-Ts(:)/647.096).^3.5 - 15.9618719*(1-Ts(:)/647.096).^4 ...
      + 1.80122502*(1-Ts(:)/647.096).^7.5));
e = max(abs(sa.P ./ aux - 1));
[np,nf] = check(sprintf('saturation vs IAPWS-95 aux eq: max %.1e', e), e < 2e-4, np, nf);
[np,nf] = check('saturation at 373.124 K ~ 0.101325 MPa', abs(sa.P(3) - 0.101325) < 2e-5, np, nf);
Tb = [230 240 250 260 270 273.16];
su = SF_coexistence('sublimation', Tb);
th = Tb(:) / 273.16;
r1408 = 611.657e-6 * exp((-0.212144006e2*th.^0.333333333e-2 + 0.273203819e2*th.^0.120666667e1 ...
        - 0.610598130e1*th.^0.170333333e1) ./ th);
e = max(abs(su.P ./ r1408 - 1));
[np,nf] = check(sprintf('sublimation (Ih + water_Brown2026) vs IAPWS R14-08: max %.1e', e), e < 2e-4, np, nf);
ex = SF_coexistence('sublimation', [175 200]);                          % default: extension on
nx = SF_coexistence('sublimation', [175 200], 'dilute_extension', false);
th = [175; 200] / 273.16;
r1408 = 611.657e-6 * exp((-0.212144006e2*th.^0.333333333e-2 + 0.273203819e2*th.^0.120666667e1 ...
        - 0.610598130e1*th.^0.170333333e1) ./ th);
[np,nf] = check('sublimation below 230 K: dilute extension by default (R14-08 +-2e-4), NaN when off', ...
    all(isnan(nx.P)) && max(abs(ex.P ./ r1408 - 1)) < 2e-4, np, nf);

% ---- P -> rho: the stable branch only returns thermodynamically stable
%      roots (Cv > 0, dP/drho > 0), also near the small (dP/drho)_T loops
%      inside the dome near Tc ------------------------------------------------
gu = SF_getprop({linspace(5, 120, 47), linspace(600, 660, 61)}, 'water_Brown2026', {'rho','Cv','Kt'});
[np,nf] = check('near-critical grid: every stable root has Cv > 0 and Kt > 0', ...
    all(isfinite(gu.rho(:))) && all(gu.Cv(:) > 0) && all(gu.Kt(:) > 0), np, nf);

fprintf('\n%d passed, %d failed\n', np, nf);
if nf > 0, error('test_fnFval:failed', '%d test(s) failed', nf); end
end

function [np, nf] = check(name, cond, np, nf)
if cond
    fprintf('[PASS] %s\n', name); np = np + 1;
else
    fprintf('[FAIL] %s\n', name); nf = nf + 1;
end
end
