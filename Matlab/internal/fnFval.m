function Results = fnFval(sp, input, props, mode, branch)
% FNFVAL  Evaluate thermodynamic properties from a Helmholtz energy spline.
%
%   Results = fnFval(sp, PT)
%   Results = fnFval(sp, PT, props)
%   Results = fnFval(sp, rhoT, props, 'rhoT')
%   Results = fnFval(sp, PT, props, 'PT', branch)
%
%   Counterpart of fnGval for equations of state given as a Helmholtz energy
%   F(rho,T) (J/kg) in density (kg/m^3) and temperature (K).  Output fields,
%   names and units are identical to fnGval so SF_getprop and its callers can
%   use either evaluator interchangeably.
%
%   Two representations are supported, selected by sp.eos:
%     'F_rhoT'  a plain tensor B-spline of F with knots{1} = rho, knots{2} = T
%               (optional sp.Tc: spline in tau = log(T/Tc) instead of T)
%     'psi'     the psi-spline surface of lbf-thermo (residual dimensionless
%               Helmholtz energy plus reference terms), evaluated by psi_val
%
%   INPUT
%     sp    - spline struct (see above)
%     input - evaluation points:
%               mode 'PT'   (default): {P,T} grid or [P(:) T(:)] scatter,
%                                      P in MPa.  Density is found by
%                                      solving P = rho^2 dF/drho on a
%                                      mechanically stable branch (below).
%               mode 'rhoT'          : {rho,T} grid or [rho(:) T(:)].
%     props - cell array of property names (or a single string); omit or
%             pass [] for all.  Same names as fnGval's pure-phase set.
%     mode  - 'PT' or 'rhoT'.
%     branch- PT mode only, which root of P(rho,T) = P to return when the
%             fluid has both a vapour-like and a liquid-like stable root:
%               'stable' (default) the one with the lower Gibbs energy
%                        (vapour below the saturation pressure, liquid above)
%               'liquid' the densest root (metastable superheated liquid
%                        below the saturation pressure)
%               'vapor'  the least dense root (metastable supersaturated
%                        vapour above the saturation pressure)
%
%   OUTPUT
%     Results - struct with the requested fields.  Grids are n1-by-nT,
%               scatter inputs return columns.  P (MPa) and rho are both
%               available in either mode: the input coordinate is echoed
%               (as a column for grids, like fnGval), the other is computed.
%               NaN where the point is outside the spline domain or, in PT
%               mode, where no mechanically stable density solves P(rho,T) = P.
%
%   THERMODYNAMICS (SI inside, P and moduli returned in MPa)
%     P     = rho^2 F_r               G  = F + rho F_r       A = F
%     S     = -F_T                    U  = F + T S           H = U + P/rho
%     Cv    = -T F_TT                 Kt = rho dP/drho
%     alpha = (dP/dT)_rho / Kt        Cp = Cv + T (dP/dT)^2 / (rho^2 dP/drho)
%     Ks    = Kt Cp/Cv                vel = sqrt(Ks/rho)
%     Kp    = 1 + rho (d2P/drho2) / (dP/drho)
%     with dP/drho = 2 rho F_r + rho^2 F_rr,   dP/dT = rho^2 F_rT,
%          d2P/drho2 = 2 F_r + 4 rho F_rr + rho^2 F_rrr.
%
%   See also: fnGval, psi_val, sp_val, SF_getprop.

    if nargin < 3, props = []; end
    if nargin < 4 || isempty(mode), mode = 'PT'; end
    if nargin < 5 || isempty(branch), branch = 'stable'; end
    mode = char(mode); branch = char(branch);
    if ~any(strcmp(branch, {'stable', 'liquid', 'vapor'}))
        error('fnFval:badInput', 'branch must be ''stable'', ''liquid'' or ''vapor'' (got ''%s'').', branch);
    end
    if ~any(strcmp(mode, {'PT', 'rhoT'}))
        error('fnFval:badInput', 'mode must be ''PT'' or ''rhoT'' (got ''%s'').', mode);
    end
    if numel(sp.knots) ~= 2
        error('fnFval:badSpline', 'Helmholtz splines must be 2-D F(rho,T).');
    end

    MPa2Pa = 1e6;

    % ------------------------------------------------------------------
    % Requested properties
    % ------------------------------------------------------------------
    defs = sf_material_defs();
    all_props = defs.base_props;
    if isempty(props)
        req = all_props;
    else
        if ischar(props) || isstring(props), props = cellstr(props); end
        unknown = setdiff(props, all_props);
        if ~isempty(unknown)
            error('fnFval:unknownProperty', ...
                  'fnFval: unsupported property name(s): %s\nValid names: %s', ...
                  strjoin(unknown, ', '), strjoin(all_props, ', '));
        end
        req = props;
    end
    want = struct();
    for ii = 1:numel(all_props), want.(all_props{ii}) = false; end
    for ii = 1:numel(req),       want.(req{ii})       = true;  end

    % Properties built on Cp / alpha share the (dP/dT, dP/drho) pieces
    need_Cv    = want.Cv || want.Cp || want.Ks || want.vel || want.Js || want.gamma_Gruneisen;
    need_alpha = want.alpha || want.Cp || want.Ks || want.vel || want.Js || want.gamma_Gruneisen;

    need.F    = want.G || want.U || want.H || want.A;
    need.Fr   = true;                              % P itself; always cheap
    need.Frr  = need_alpha || want.Kt || want.Kp;
    need.Frrr = want.Kp;
    need.FT   = want.S || want.U || want.H;
    need.FTT  = need_Cv;
    need.FrT  = need_alpha;

    % ------------------------------------------------------------------
    % Parse input
    % ------------------------------------------------------------------
    gridded = iscell(input);
    if gridded
        if numel(input) ~= 2
            error('fnFval:badInput', 'Gridded input must be a 2-element cell {X,T}.');
        end
        X = input{1}(:);  T = input{2}(:);
        [Xm, Tm] = ndgrid(X, T);
    else
        if size(input, 2) ~= 2
            error('fnFval:badInput', 'Scattered input must be N-by-2 [X T].');
        end
        Xm = input(:,1);  Tm = input(:,2);
    end
    sz = size(Xm);
    [rlim, Tlim] = helm_domain(sp);

    % ------------------------------------------------------------------
    % Density at every point
    % ------------------------------------------------------------------
    if strcmp(mode, 'PT')
        rhom = invert_P(sp, Xm(:), Tm(:), rlim, Tlim, branch);
        rhom = reshape(rhom, sz);
    else
        rhom = Xm;
        rhom(rhom < rlim(1) | rhom > rlim(2)) = NaN;
    end
    rhom(Tm < Tlim(1) | Tm > Tlim(2)) = NaN;

    % ------------------------------------------------------------------
    % F derivatives at the valid (rho,T) points
    % ------------------------------------------------------------------
    ok = isfinite(rhom);
    d  = helm_derivs(sp, rhom(ok), Tm(ok), need);
    names = fieldnames(d);
    for ii = 1:numel(names)
        v = NaN(sz);
        v(ok) = d.(names{ii});
        d.(names{ii}) = v;
    end

    % ------------------------------------------------------------------
    % Properties
    % ------------------------------------------------------------------
    r  = rhom;
    Pp = r.^2 .* d.Fr;                              % Pa
    if need.Frr
        dPdr = 2*r.*d.Fr + r.^2.*d.Frr;             % Pa m^3/kg
    end
    if need_Cv,    Cv    = -Tm .* d.FTT;            end
    if need_alpha
        dPdT  = r.^2 .* d.FrT;                      % Pa/K
        alpha = dPdT ./ (r .* dPdr);
        Cp    = [];
        if need_Cv, Cp = Cv + Tm .* dPdT.^2 ./ (r.^2 .* dPdr); end
    end
    Kt_Pa = [];
    if need.Frr, Kt_Pa = r .* dPdr; end

    Results = struct();
    if want.G,     Results.G     = d.F + Pp ./ r;                      end
    if want.S,     Results.S     = -d.FT;                              end
    if want.U,     Results.U     = d.F - Tm .* d.FT;                   end
    if want.H,     Results.H     = d.F - Tm .* d.FT + Pp ./ r;         end
    if want.A,     Results.A     = d.F;                                end
    if want.rho
        if strcmp(mode, 'rhoT') && gridded
            Results.rho = input{1}(:);
        else
            Results.rho = rhom;
        end
    end
    if want.Cp,    Results.Cp    = Cp;                                 end
    if want.Cv,    Results.Cv    = Cv;                                 end
    if want.Kt,    Results.Kt    = Kt_Pa / MPa2Pa;                     end
    if want.Ks,    Results.Ks    = Kt_Pa .* Cp ./ Cv / MPa2Pa;         end
    if want.Kp
        d2Pdr2 = 2*d.Fr + 4*r.*d.Frr + r.^2.*d.Frrr;
        Results.Kp = 1 + r .* d2Pdr2 ./ dPdr;
    end
    if want.alpha, Results.alpha = alpha;                              end
    if want.vel
        w2 = Kt_Pa .* Cp ./ Cv ./ r;
        w2(w2 < 0) = NaN;                           % mechanically unstable
        Results.vel = sqrt(w2);
    end
    if want.Js,    Results.Js    = Tm .* alpha ./ (r .* Cp) * MPa2Pa;  end
    if want.gamma_Gruneisen
        Results.gamma_Gruneisen = alpha .* Kt_Pa ./ (r .* Cv);
    end
    if want.P
        if strcmp(mode, 'PT')
            if gridded, Results.P = input{1}(:); else, Results.P = Xm; end
        else
            Results.P = Pp / MPa2Pa;
        end
    end
    if want.T
        if gridded, Results.T = input{2}(:); else, Results.T = Tm; end
    end
end


% ======================================================================
%  SUBFUNCTIONS
% ======================================================================

function tf = is_psi(sp)
    tf = isfield(sp, 'eos') && strcmp(sp.eos, 'psi');
end


function [rlim, Tlim] = helm_domain(sp)
% HELM_DOMAIN  Valid density and temperature range of the spline.
    if is_psi(sp)
        % below the lowest density knot the psi surface is continued (virial)
        rlim = [0, sp.rhoc * exp(3 * sp.knots{1}(end))];
        Tlim = sp.Tc * exp([sp.knots{2}(1), sp.knots{2}(end)]);
    else
        rlim = [sp.knots{1}(1), sp.knots{1}(end)];
        if isfield(sp, 'Tc')
            Tlim = sp.Tc * exp([sp.knots{2}(1), sp.knots{2}(end)]);
        else
            Tlim = [sp.knots{2}(1), sp.knots{2}(end)];
        end
    end
end


function d = helm_derivs(sp, r, T, need)
% HELM_DERIVS  F and its (rho,T) derivatives at scattered points.
    fl = {'F','Fr','Frr','Frrr','FT','FTT','FrT'};
    for k = 1:numel(fl)
        if ~isfield(need, fl{k}), need.(fl{k}) = false; end
    end
    if is_psi(sp)
        d = psi_val(sp, r, T, need);
    else
        d = spline_derivs(sp, r, T, need);
    end
end


function d = spline_derivs(sp, r, T, need)
% SPLINE_DERIVS  Plain F(rho,T) B-spline, derivatives in T units.
    flgTc = isfield(sp, 'Tc');
    if flgTc, x = [r(:), log(T(:) / sp.Tc)]; else, x = [r(:), T(:)]; end
    n = numel(r);
    d = struct();
    if n == 0
        fl = {'F','Fr','Frr','Frrr','FT','FTT','FrT'};
        for k = 1:numel(fl), if need.(fl{k}), d.(fl{k}) = zeros(0, 1); end, end
        return
    end

    if need.F,    d.F    = sp_val(sp, [0 0], x); end
    if need.Fr,   d.Fr   = sp_val(sp, [1 0], x); end
    if need.Frr,  d.Frr  = sp_val(sp, [2 0], x); end
    if need.Frrr, d.Frrr = sp_val(sp, [3 0], x); end
    need_FT = need.FT || (flgTc && need.FTT);
    if need_FT,   d.FT   = sp_val(sp, [0 1], x); end
    if need.FTT,  d.FTT  = sp_val(sp, [0 2], x); end
    if need.FrT,  d.FrT  = sp_val(sp, [1 1], x); end

    if flgTc
        % d/dT = (1/T) d/dtau ;  d2/dT2 = (d2/dtau2 - d/dtau) / T^2
        Tc = T(:);
        if need.FTT, d.FTT = (d.FTT - d.FT) ./ Tc.^2; end
        if need_FT,  d.FT  = d.FT ./ Tc;              end
        if need.FrT, d.FrT = d.FrT ./ Tc;             end
    end
    if need_FT && ~need.FT, d = rmfield(d, 'FT'); end
end


function [P, dPdr] = P_dPdrho(sp, r, T)
% P (Pa) and (dP/drho)_T (Pa m^3/kg) at scattered points.
    nd = struct('Fr', true, 'Frr', true);
    d = helm_derivs(sp, r(:), T(:), nd);
    P    = r(:).^2 .* d.Fr;
    dPdr = 2 * r(:) .* d.Fr + r(:).^2 .* d.Frr;
end


function rho = invert_P(sp, P, T, rlim, Tlim, branch)
% INVERT_P  Density solving rho^2 dF/drho = P at each point (P in MPa).
%
%   Brackets every root on a coarse density grid (per distinct temperature)
%   where P increases with rho (mechanically stable), refines the least and
%   the most dense of them with a bracket-safeguarded Newton iteration, and
%   returns the one selected by branch ('stable': lower Gibbs energy).

    n   = numel(P);
    rho = NaN(n, 1);
    inT = isfinite(P) & isfinite(T) & T >= Tlim(1) & T <= Tlim(2);
    if ~any(inT), return; end

    Ppa = P * 1e6;

    % ---- coarse brackets, one pass per distinct T -----------------------
    nr = 600;
    if is_psi(sp)
        % geometric spacing (the surface spans decades), extended below the
        % lowest density knot into the virial continuation (dilute vapour)
        x_lo = min(sp.knots{1}(1), log(1e-12 / sp.rhoc) / 3);
        rg = sp.rhoc * exp(3 * linspace(x_lo, sp.knots{1}(end), nr)).';
        % plus a fine uniform sampling of the dense fluid, so that the small
        % (dP/drho)_T loops inside the dome near Tc never share a bracket
        % with the true liquid or vapour root
        rtop = sp.rhoc * exp(3 * sp.knots{1}(end));
        rg = unique([rg; (1:2:min(2000, rtop)).']);
        nr = numel(rg);
    else
        rg = linspace(rlim(1), rlim(2), nr).';
    end
    klo = zeros(n, 1);  khi = zeros(n, 1);
    idx  = find(inT);
    [Tu, ~, jT] = unique(T(idx));
    nTu = numel(Tu);
    if is_psi(sp)
        dg = psi_val(sp, rg, Tu, struct('Fr', true, 'Frr', true), 'grid');   % tensor grid: fast
        Fr = dg.Fr; Frr = dg.Frr;
    else
        if isfield(sp, 'Tc'), tu = log(Tu / sp.Tc); else, tu = Tu; end
        Fr  = reshape(sp_val(sp, [1 0], {rg.', tu(:).'}), nr, nTu);
        Frr = reshape(sp_val(sp, [2 0], {rg.', tu(:).'}), nr, nTu);
    end
    pg  = bsxfun(@times, rg.^2, Fr);                            % P, Pa
    dpg = bsxfun(@times, 2 * rg, Fr) + bsxfun(@times, rg.^2, Frr);   % dP/drho
    col = zeros(n, 1); col(idx) = jT;                            % T column of each point
    for j = 1:nTu
        pts = idx(jT == j);
        s   = pg(:, j) - Ppa(pts).';                 % nr-by-npts
        up  = s(1:end-1, :) <= 0 & s(2:end, :) >= 0 & diff(pg(:, j)) > 0;
        kk  = bsxfun(@times, up, (1:nr-1).');
        khi(pts) = max(kk, [], 1);                   % densest crossing
        kk(~up) = Inf;
        kl  = min(kk, [], 1); kl(~isfinite(kl)) = 0;
        klo(pts) = kl;                               % least dense crossing
    end
    if strcmp(branch, 'liquid'), klo = khi; end
    if strcmp(branch, 'vapor'),  khi = klo; end

    % ---- refine the candidate roots ----------------------------------------
    hasr = khi > 0;
    ih = find(hasr);
    rho_hi = NaN(n, 1);
    x0 = hermite_guess(rg, pg, dpg, khi(ih), col(ih), Ppa(ih));
    rho_hi(ih) = newton_bracket(sp, rg(khi(ih)), rg(khi(ih) + 1), Ppa(ih), T(ih), x0);
    rho = rho_hi;
    two = find(hasr & klo ~= khi);
    if ~isempty(two)
        x0 = hermite_guess(rg, pg, dpg, klo(two), col(two), Ppa(two));
        rlo = newton_bracket(sp, rg(klo(two)), rg(klo(two) + 1), Ppa(two), T(two), x0);
        % lower Gibbs energy wins: G = F + rho F_r
        nd = struct('F', true, 'Fr', true);
        dh = helm_derivs(sp, rho_hi(two), T(two), nd);
        dl = helm_derivs(sp, rlo, T(two), nd);
        Gh = dh.F + rho_hi(two) .* dh.Fr;
        Gl = dl.F + rlo .* dl.Fr;
        usev = Gl < Gh | ~isfinite(Gh);
        rho(two(usev)) = rlo(usev);
    end
end


function x = hermite_guess(rg, pg, dpg, k, j, Pt)
% Starting density inside each bracket [rg(k), rg(k+1)] at T column j: the
% root of the cubic Hermite interpolant of P(rho) built from P and dP/drho at
% the bracket ends (already known from the grid), so Newton needs ~1 step.
    k = k(:); j = j(:); Pt = Pt(:);
    nr = numel(rg);
    i0 = k + (j - 1) * nr; i1 = i0 + 1;
    a = rg(k); h = rg(k + 1) - a;
    pa = pg(i0); pb = pg(i1); da = dpg(i0) .* h; db = dpg(i1) .* h;
    t = min(max((Pt - pa) ./ (pb - pa), 0), 1);
    t(~isfinite(t)) = 0.5;
    for it = 1:12
        t2 = t.^2; t3 = t2 .* t;
        f  = (2*t3 - 3*t2 + 1) .* pa + (t3 - 2*t2 + t) .* da + (-2*t3 + 3*t2) .* pb + (t3 - t2) .* db - Pt;
        df = (6*t2 - 6*t) .* pa + (3*t2 - 4*t + 1) .* da + (-6*t2 + 6*t) .* pb + (3*t2 - 2*t) .* db;
        tn = t - f ./ df;
        bad = ~isfinite(tn) | tn < 0 | tn > 1;
        tn(bad) = t(bad);
        t = tn;
    end
    x = a + t .* h;
end


function x = newton_bracket(sp, a, b, Pt, ta, x0)
% Bracket-safeguarded Newton for rho^2 F_r = Pt (Pa) with P increasing on [a,b].
% Converged when the Newton step is below 1e-9 rho: with quadratic
% convergence the density after that step is accurate to ~1e-16 relative.
    a = a(:); b = b(:); Pt = Pt(:); ta = ta(:);
    if nargin >= 6 && ~isempty(x0)
        x = x0(:);
    else
        pa = P_dPdrho(sp, a, ta);
        pb = P_dPdrho(sp, b, ta);
        x  = a + (Pt - pa) .* (b - a) ./ (pb - pa);
    end
    bad = ~isfinite(x) | x < a | x > b;
    x(bad) = (a(bad) + b(bad)) / 2;
    live = true(numel(x), 1);
    for it = 1:80
        [f, df] = P_dPdrho(sp, x(live), ta(live));
        f = f - Pt(live);
        al = a(live); bl = b(live); xl = x(live);
        % shrink bracket with the sign of the residual (P increasing in rho)
        al(f < 0) = xl(f < 0);
        bl(f > 0) = xl(f > 0);
        xn  = xl - f ./ df;
        out = ~isfinite(xn) | xn <= al | xn >= bl;
        xn(out) = (al(out) + bl(out)) / 2;
        step = abs(xn - xl);
        a(live) = al; b(live) = bl; x(live) = xn;
        done = (step <= 1e-9 * xn & ~out) | f == 0;
        live(live) = ~done;
        if ~any(live), break; end
    end
end
