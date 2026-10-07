function d = psi_val(sp, rho, T, need, mode)
% PSI_VAL  Helmholtz energy and its (rho,T) derivatives from a psi-spline surface.
%
%   d = psi_val(sp, rho, T, need)
%   d = psi_val(sp, rho, T, need, 'grid')   % tensor grid rho x T, F_r / F_rr only
%
%   Toolbox-free port of psiH2O_val.m (lbf-thermo, JMB 2026) for the
%   "psi" representation of pure water: the residual dimensionless Helmholtz
%   energy is a tensor B-spline in stretched coordinates plus reference terms,
%
%     phi(delta,tau) = phi0 + dphi_ref + phi_crit + phi_2s + delta * psi(x,y)
%     x = ln(rho/rhoc)/3,  y = ln(T/Tc),  delta = rho/rhoc,  tau = Tc/T,
%     F = R T phi   (J/kg)
%
%   with
%     phi0      Planck-Einstein ideal gas (sp.phi0_n0, sp.phi0_g0)
%     dphi_ref  reacting-mixture reference table, B-spline in (x,y) clamped at
%               its edges, with a low-density extension sp.dphi_ref.low joined
%               C^1 below the table floor
%     phi_crit  tabulated KW2000 critical term sp.dphi_ref.crit_table, a
%               B-spline in (1 - tau, delta - 1), zero outside its box
%     phi_2s    analytic low-T two-structure term (sp.dphi_ref.lowT2s_params)
%     psi       sp itself; below its lowest x knot it is continued as a virial
%               form (C^1 at the floor)
%
%   INPUT
%     sp    - the 'sp_psi' struct of a psi surface (SeaFreeze stores it as sp
%             with sp.eos = 'psi').  Required: knots, coefs, number, order, dim,
%             Tc, rhoc, R.  Optional: phi0_n0/phi0_g0, dphi_ref (with low,
%             crit_table, lowT2s_params).
%     rho,T - column vectors of scattered states (kg/m^3, K)
%     need  - struct of logicals: F, Fr, Frr, Frrr, FT, FTT, FrT (absent = false)
%     mode  - 'grid': rho and T are the axes of a tensor grid; outputs are
%             numel(rho)-by-numel(T).  Every spline is then evaluated with
%             sp_val's gridded (per-dimension de Boor) path, which is far
%             cheaper than scattered evaluation.  Only Fr and Frr are
%             supported (used to bracket the roots of P(rho,T) = P).
%
%   OUTPUT
%     d - struct with the requested fields, columns like rho, NaN where the
%         state lies outside the surface (x above the top density knot, or T
%         outside the T knots).  Units: F J/kg; Fr J m^3/kg^2; Frr J m^6/kg^3;
%         FT J/(kg K); FTT J/(kg K^2); FrT J m^3/(kg^2 K).
%         Frrr is a central finite difference of Frr in rho (the reference
%         terms carry no analytic third derivative).
%
%   Mirrors psiH2O_val.m field by field (fnval/fnder -> sp_val).  The retired
%   Stage 5-6 two-state terms (twostate_params) are not ported; a surface
%   carrying them is refused.
%
%   See also: fnFval, sp_val.

    fl = {'F','Fr','Frr','Frrr','FT','FTT','FrT'};
    for k = 1:numel(fl)
        if ~isfield(need, fl{k}), need.(fl{k}) = false; end
    end
    if nargin >= 5 && strcmp(mode, 'grid')
        if need.F || need.Frrr || need.FT || need.FTT || need.FrT
            error('psi_val:grid', 'grid mode supports Fr and Frr only.');
        end
        d = psi_grid_rho(sp, rho(:), T(:), need);
        return
    end

    rho = rho(:); T = T(:); n = numel(rho);
    Tc = double(sp.Tc); rhoc = double(sp.rhoc); R = double(sp.R);

    has_ref = isfield(sp, 'dphi_ref');
    if (isfield(sp, 'twostate_params') || isfield(sp, 'twostate_s_params')) || ...
       (has_ref && (isfield(sp.dphi_ref, 'twostate_params') || isfield(sp.dphi_ref, 'twostate_s_params')))
        error('psi_val:unsupported', 'retired two-state terms are not supported by psi_val.');
    end
    if has_ref && isfield(sp.dphi_ref, 'crit_params') && ~isfield(sp.dphi_ref, 'crit_table')
        error('psi_val:unsupported', 'direct KW2000 critical term is not supported; the surface must carry crit_table.');
    end

    % ---- third derivative by finite difference (recursive) --------------
    if need.Frrr
        h = 1e-4;
        nd = struct('Frr', true);
        dp = psi_val(sp, rho * (1 + h), T, nd);
        dm = psi_val(sp, rho * (1 - h), T, nd);
        need.Frrr = false;
        d = psi_val(sp, rho, T, need);
        d.Frrr = (dp.Frr - dm.Frr) ./ (2 * h * rho);
        return
    end

    % ---- which phi derivatives are required -------------------------------
    wp  = need.F || need.FT;
    wd  = need.Fr || need.FrT;
    wdd = need.Frr;
    wt  = need.FT;
    wtt = need.FTT;
    wdt = need.FrT;
    wy  = wt || wtt || wdt;

    x = log(rho / rhoc) / 3;  y = log(T / Tc);
    kx = sp.knots{1}(:); ky = sp.knots{2}(:);
    ok = isfinite(x) & isfinite(y) & rho > 0 & x <= kx(end) & y >= ky(1) & y <= ky(end);

    del = rho / rhoc;  tau = Tc ./ T;
    z = NaN(n, 1);
    phir = z; phir_d = z; phir_dd = z; phir_t = z; phir_tt = z; phir_dt = z;

    if any(ok)
        xo = x(ok); yo = y(ok); dlt = del(ok); to = tau(ok);
        pts = [xo, yo];

        % ---- spline residual psi and derivatives in (x,y) -----------------
        s   = sp_val(sp, [0 0], pts);
        sx  = sp_val(sp, [1 0], pts);
        sxx = []; sy = []; syy = []; sxy = [];
        if wdd, sxx = sp_val(sp, [2 0], pts); end
        if wy,  sy  = sp_val(sp, [0 1], pts); end
        if wtt, syy = sp_val(sp, [0 2], pts); end
        if wdt, sxy = sp_val(sp, [1 1], pts); end

        % virial continuation below the lowest density knot (C^1 join)
        x0 = kx(1);
        below = xo < x0;
        if any(below)
            qb = [repmat(x0, nnz(below), 1), yo(below)];
            r  = exp(3 * (xo(below) - x0));
            s00 = sp_val(sp, [0 0], qb); s10 = sp_val(sp, [1 0], qb);
            s(below)  = s00 + s10 .* (r - 1) / 3;
            sx(below) = s10 .* r;
            if wdd, sxx(below) = 3 * s10 .* r; end
            if wy
                s01 = sp_val(sp, [0 1], qb); s11 = sp_val(sp, [1 1], qb);
                sy(below) = s01 + s11 .* (r - 1) / 3;
                if wtt
                    s02 = sp_val(sp, [0 2], qb); s12 = sp_val(sp, [1 2], qb);
                    syy(below) = s02 + s12 .* (r - 1) / 3;
                end
                if wdt, sxy(below) = s11 .* r; end
            end
        end

        pr    = dlt .* s;
        pr_d  = s + sx / 3;
        pr_dd = []; pr_t = []; pr_tt = []; pr_dt = [];
        if wdd, pr_dd = (sx / 3 + sxx / 9) ./ dlt;        end
        if wt,  pr_t  = -dlt .* sy ./ to;                  end
        if wtt, pr_tt =  dlt .* (syy + sy) ./ to.^2;       end
        if wdt, pr_dt = -(sy + sxy / 3) ./ to;            end

        % ---- reference table dphi_ref (clamped; low-density extension) ---
        if has_ref
            rs = sp.dphi_ref;
            rkx = rs.knots{1}(:); rky = rs.knots{2}(:);
            Xc = min(max(xo, rkx(1)), rkx(end)); Yc = min(max(yo, rky(1)), rky(end));
            oob = (xo ~= Xc) | (yo ~= Yc);
            pr2 = [Xc, Yc];
            f = sp_val(rs, [0 0], pr2);
            fx = sp_val(rs, [1 0], pr2); fx(oob) = 0;
            fxx = []; fy = []; fyy = []; fxy = [];
            if wdd, fxx = sp_val(rs, [2 0], pr2); fxx(oob) = 0; end
            if wy,  fy  = sp_val(rs, [0 1], pr2); fy(oob)  = 0; end
            if wtt, fyy = sp_val(rs, [0 2], pr2); fyy(oob) = 0; end
            if wdt, fxy = sp_val(rs, [1 1], pr2); fxy(oob) = 0; end
            if isfield(rs, 'low')
                lw = rs.low; lkx = lw.knots{1}(:);
                bl = xo < rkx(1) & xo >= lkx(1) & yo >= rky(1) & yo <= rky(end);
                if any(bl)
                    ql = [xo(bl), yo(bl)]; q0 = [repmat(rkx(1), nnz(bl), 1), yo(bl)];
                    r  = exp(3 * (xo(bl) - rkx(1)));
                    L = @(i, j) sp_val(lw, [i j], ql);
                    D = @(i, j) sp_val(rs, [i j], q0) - sp_val(lw, [i j], q0);
                    D10 = D(1, 0);
                    f(bl)  = L(0, 0) + D(0, 0) + D10 .* (r - 1) / 3;
                    fx(bl) = L(1, 0) + D10 .* r;
                    if wdd, fxx(bl) = L(2, 0) + 3 * D10 .* r; end
                    if wy
                        D01 = D(0, 1); D11 = D(1, 1);
                        fy(bl) = L(0, 1) + D01 + D11 .* (r - 1) / 3;
                        if wtt, fyy(bl) = L(0, 2) + D(0, 2) + D(1, 2) .* (r - 1) / 3; end
                        if wdt, fxy(bl) = L(1, 1) + D11 .* r; end
                    end
                end
            end
            pr   = pr   + f;
            pr_d = pr_d + fx ./ (3 * dlt);
            if wdd, pr_dd = pr_dd + (fxx / 9 - fx / 3) ./ dlt.^2; end
            if wt,  pr_t  = pr_t  - fy ./ to;                     end
            if wtt, pr_tt = pr_tt + (fyy + fy) ./ to.^2;          end
            if wdt, pr_dt = pr_dt - fxy ./ (3 * dlt .* to);        end
        end

        % ---- tabulated KW2000 critical term ---------------------------------
        ct = [];
        if isfield(sp, 'crit_table'), ct = sp.crit_table;
        elseif has_ref && isfield(sp.dphi_ref, 'crit_table'), ct = sp.dphi_ref.crit_table; end
        if ~isempty(ct)
            c0 = double(ct.c0);
            ckx = ct.knots{1}(:); cky = ct.knots{2}(:);
            xq = 1 - to; yq = dlt - 1;
            inb = xq >= ckx(1) & xq <= ckx(end) & yq >= cky(1) & yq <= cky(end);
            if any(inb)
                pq = [xq(inb), yq(inb)];
                g = zeros(nnz(ok), 1); gy = g; gyy = g; gx = g; gxx = g; gxy = g;
                g(inb)  = sp_val(ct, [0 0], pq);
                gy(inb) = sp_val(ct, [0 1], pq);
                if wdd, gyy(inb) = sp_val(ct, [0 2], pq); end
                if wy,  gx(inb)  = sp_val(ct, [1 0], pq); end
                if wtt, gxx(inb) = sp_val(ct, [2 0], pq); end
                if wdt, gxy(inb) = sp_val(ct, [1 1], pq); end
                pr   = pr   + c0 * g ./ dlt;
                pr_d = pr_d + c0 * (gy ./ dlt - g ./ dlt.^2);
                if wdd, pr_dd = pr_dd + c0 * (gyy ./ dlt - 2 * gy ./ dlt.^2 + 2 * g ./ dlt.^3); end
                if wt,  pr_t  = pr_t  - c0 * gx ./ dlt;                    end
                if wtt, pr_tt = pr_tt + c0 * gxx ./ dlt;                   end
                if wdt, pr_dt = pr_dt - c0 * (gxy ./ dlt - gx ./ dlt.^2);   end
            end
        end

        % ---- low-T two-structure term (analytic) ----------------------------
        lj = '';
        if isfield(sp, 'lowT2s_params'), lj = sp.lowT2s_params;
        elseif has_ref && isfield(sp.dphi_ref, 'lowT2s_params'), lj = sp.dphi_ref.lowT2s_params; end
        if ~isempty(lj)
            q = lowT2s_phi(dlt, to, lowT2s_decode(char(lj)), R);
            pr   = pr   + q.p;
            pr_d = pr_d + q.d;
            if wdd, pr_dd = pr_dd + q.dd; end
            if wt,  pr_t  = pr_t  + q.t;  end
            if wtt, pr_tt = pr_tt + q.tt; end
            if wdt, pr_dt = pr_dt + q.dt; end
        end

        phir(ok) = pr; phir_d(ok) = pr_d;
        if wdd, phir_dd(ok) = pr_dd; end
        if wt,  phir_t(ok)  = pr_t;  end
        if wtt, phir_tt(ok) = pr_tt; end
        if wdt, phir_dt(ok) = pr_dt; end
    end

    % ---- ideal-gas part ------------------------------------------------------
    if isfield(sp, 'phi0_n0') && isfield(sp, 'phi0_g0')
        n0 = double(sp.phi0_n0(:)'); g0 = double(sp.phi0_g0(:)');
    else
        n0 = [-8.3204464837497 6.6832105275932 3.00632 0.012436 0.97315 1.27950 0.96956 0.24873];
        g0 = [1.28728967 3.53734222 7.74073708 9.24437796 27.5075105];
    end
    phi0 = log(del) + n0(1) + n0(2) * tau + n0(3) * log(tau);
    phi0_t  = n0(2) + n0(3) ./ tau;
    phi0_tt = -n0(3) ./ tau.^2;
    for i = 1:numel(g0)
        e = exp(-g0(i) * tau);
        phi0    = phi0    + n0(3+i) * log(1 - e);
        phi0_t  = phi0_t  + n0(3+i) * g0(i) * (1 ./ (1 - e) - 1);
        phi0_tt = phi0_tt - n0(3+i) * g0(i)^2 * e ./ (1 - e).^2;
    end

    % ---- F and its (rho,T) derivatives -----------------------------------------
    %   F = R T phi;  d/drho = (1/rhoc) d/ddelta;  d tau/dT = -tau/T
    RT = R * T;
    d = struct();
    if wp
        phi = phi0 + phir;
        if need.F,  d.F  = RT .* phi;                                  end
        if need.FT, d.FT = R * (phi - tau .* (phi0_t + phir_t));       end
    end
    if wd
        phi_d = 1 ./ del + phir_d;
        if need.Fr,  d.Fr  = RT .* phi_d / rhoc;                        end
        if need.FrT, d.FrT = R / rhoc * (phi_d - tau .* phir_dt);        end
    end
    if need.Frr, d.Frr = RT .* (-1 ./ del.^2 + phir_dd) / rhoc^2;         end
    if need.FTT, d.FTT = R * tau.^2 .* (phi0_tt + phir_tt) ./ T;          end
    fn = fieldnames(d);
    for k = 1:numel(fn), v = d.(fn{k}); v(~ok) = NaN; d.(fn{k}) = v; end
end


% ==========================================================================
function d = psi_grid_rho(sp, rv, Tv, need)
% F_r and F_rr on the tensor grid rv (nr) x Tv (nT); the same terms as the
% scattered path, with every spline evaluated on its grid.
    Tc = double(sp.Tc); rhoc = double(sp.rhoc); R = double(sp.R);
    nr = numel(rv); nT = numel(Tv);
    EV = @(spl, dv, xs, ys) reshape(sp_val(spl, dv, {xs(:).', ys(:).'}), numel(xs), numel(ys));
    x = log(rv / rhoc) / 3;                 % nr x 1
    y = log(Tv(:).' / Tc);                  % 1 x nT
    DEL = repmat(rv / rhoc, 1, nT);
    TAU = repmat(Tc ./ Tv(:).', nr, 1);
    kx = sp.knots{1}(:); ky = sp.knots{2}(:);
    ok = (x <= kx(end)) & (y >= ky(1) & y <= ky(end));      % nr x nT

    % spline residual, with the virial continuation below the lowest x knot
    s = EV(sp, [0 0], x, y); sx = EV(sp, [1 0], x, y); sxx = EV(sp, [2 0], x, y);
    below = x < kx(1);
    if any(below)
        s0 = EV(sp, [0 0], kx(1), y); s1 = EV(sp, [1 0], kx(1), y);
        r = exp(3 * (x(below) - kx(1)));
        s(below, :)   = bsxfun(@plus, s0, bsxfun(@times, s1, (r - 1) / 3));
        sx(below, :)  = bsxfun(@times, s1, r);
        sxx(below, :) = bsxfun(@times, 3 * s1, r);
    end
    pr_d  = s + sx / 3;
    pr_dd = (sx / 3 + sxx / 9) ./ DEL;

    % reference table (clamped) + low-density extension
    if isfield(sp, 'dphi_ref')
        rs = sp.dphi_ref; rkx = rs.knots{1}(:); rky = rs.knots{2}(:);
        Xc = min(max(x, rkx(1)), rkx(end)); Yc = min(max(y, rky(1)), rky(end));
        oob = bsxfun(@or, x ~= Xc, y ~= Yc);
        fx = EV(rs, [1 0], Xc, Yc); fxx = EV(rs, [2 0], Xc, Yc);
        fx(oob) = 0; fxx(oob) = 0;
        if isfield(rs, 'low')
            lw = rs.low; lkx = lw.knots{1}(:);
            blr = x < rkx(1) & x >= lkx(1);
            blc = y >= rky(1) & y <= rky(end);
            if any(blr) && any(blc)
                r = exp(3 * (x(blr) - rkx(1)));
                D10 = EV(rs, [1 0], rkx(1), y(blc)) - EV(lw, [1 0], rkx(1), y(blc));
                fx(blr, blc)  = EV(lw, [1 0], x(blr), y(blc)) + bsxfun(@times, D10, r);
                fxx(blr, blc) = EV(lw, [2 0], x(blr), y(blc)) + bsxfun(@times, 3 * D10, r);
            end
        end
        pr_d  = pr_d  + fx ./ (3 * DEL);
        pr_dd = pr_dd + (fxx / 9 - fx / 3) ./ DEL.^2;
    end

    % tabulated KW2000 critical term: table axes (1 - tau) per T, (delta - 1) per rho
    ct = [];
    if isfield(sp, 'crit_table'), ct = sp.crit_table;
    elseif isfield(sp, 'dphi_ref') && isfield(sp.dphi_ref, 'crit_table'), ct = sp.dphi_ref.crit_table; end
    if ~isempty(ct)
        c0 = double(ct.c0); ckx = ct.knots{1}(:); cky = ct.knots{2}(:);
        xq = 1 - Tc ./ Tv(:).';  yq = rv / rhoc - 1;
        inc = xq >= ckx(1) & xq <= ckx(end);
        inr = yq >= cky(1) & yq <= cky(end);
        if any(inc) && any(inr)
            g = zeros(nr, nT); gy = g; gyy = g;
            g(inr, inc)   = EV(ct, [0 0], xq(inc), yq(inr)).';
            gy(inr, inc)  = EV(ct, [0 1], xq(inc), yq(inr)).';
            gyy(inr, inc) = EV(ct, [0 2], xq(inc), yq(inr)).';
            pr_d  = pr_d  + c0 * (gy ./ DEL - g ./ DEL.^2);
            pr_dd = pr_dd + c0 * (gyy ./ DEL - 2 * gy ./ DEL.^2 + 2 * g ./ DEL.^3);
        end
    end

    % low-T two-structure term (analytic, elementwise)
    lj = '';
    if isfield(sp, 'lowT2s_params'), lj = sp.lowT2s_params;
    elseif isfield(sp, 'dphi_ref') && isfield(sp.dphi_ref, 'lowT2s_params'), lj = sp.dphi_ref.lowT2s_params; end
    if ~isempty(lj)
        q = lowT2s_phi(DEL(:), TAU(:), lowT2s_decode(char(lj)), R);
        pr_d  = pr_d  + reshape(q.d,  nr, nT);
        pr_dd = pr_dd + reshape(q.dd, nr, nT);
    end

    RT = R * repmat(Tv(:).', nr, 1);
    d = struct();
    if need.Fr,  d.Fr  = RT .* (1 ./ DEL + pr_d) / rhoc;          d.Fr(~ok)  = NaN; end
    if need.Frr, d.Frr = RT .* (-1 ./ DEL.^2 + pr_dd) / rhoc^2;   d.Frr(~ok) = NaN; end
end


% ==========================================================================
function p = lowT2s_decode(txt)
% Parse the JSON parameter string of the low-T two-structure term without
% jsondecode (Octave compatibility).  Values are numbers or null.
    p = struct();
    tok = regexp(txt, '"(\w+)"\s*:\s*([^,}]+)', 'tokens');
    for k = 1:numel(tok)
        v = strtrim(tok{k}{2});
        if strcmp(v, 'null'), p.(tok{k}{1}) = [];
        else,                 p.(tok{k}{1}) = str2double(v); end
    end
end


function q = lowT2s_phi(d, tau, p, R)
% Mirror of twostructure_lowT.TwoStructureLowT.__call__ (Python) / psiH2O_val
% lowT2s_phi: phi_2s = min_x [x G + x ln x + (1-x) ln(1-x) + omega x (1-x)] and
% its delta/tau derivatives at the equilibrium x.  Module constants TC, RHOC
% are the IAPWS-95 critical point.
    TC = 647.096; RHOC = 322.0;
    if isfield(p, 'window') && ~isempty(p.window)
        error('psi_val:lowT2s', 'windowed low-T term not supported');
    end
    s_ = p.s; th = p.Tstar / TC; be = p.b / (R * TC); ds = p.rhostar / RHOC;
    la = 0; if isfield(p, 'lam') && ~isempty(p.lam), la = p.lam; end
    % log share of the lam term (lbf-thermo 2026-09-30, stage5_19eL on); 0 = pure 1/T form
    fl = 0; if isfield(p, 'lam_log') && ~isempty(p.lam_log), fl = p.lam_log; end
    L = log(d / ds);
    if isfield(p, 'rho_m') && ~isempty(p.rho_m)
        Lm = log(p.rho_m / p.rhostar); Q = L - L.^2 / (2 * Lm); Q1 = 1 - L / Lm; Q2 = -1 / Lm + 0 * L;
    else
        Q = L; Q1 = 1 + 0 * L; Q2 = 0 * L;
    end
    g   = s_ * (th * tau - 1) + la * ((1 - fl) * (1 ./ (th * tau) - 1) - fl * log(th * tau)) + be * tau .* Q;
    gd  = be * tau .* Q1 ./ d;  gdd = be * tau .* (Q2 - Q1) ./ d.^2;
    gt  = s_ * th - la * ((1 - fl) ./ (th * tau.^2) + fl ./ tau) + be * Q;
    gtt = la * ((1 - fl) * 2 ./ (th * tau.^3) + fl ./ tau.^2);  gdt = be * Q1 ./ d;
    Lw = log(d / (p.rho_top / RHOC));
    if isfield(p, 'omega_h') && ~isempty(p.omega_h) && p.omega_h > 0
        h = p.omega_h; zz = -Lw / h;
        sp_ = h * (max(zz, 0) + log1p(exp(-abs(zz)))); sg = 1 ./ (1 + exp(-zz));
        om = p.omega_top - p.omega1 * sp_; dL = p.omega1 * sg; d2L = -(p.omega1 / h) * sg .* (1 - sg);
        omd = dL ./ d; omdd = (d2L - dL) ./ d.^2;
    else
        om = p.omega_top + p.omega1 * Lw; omd = p.omega1 ./ d; omdd = -p.omega1 ./ d.^2;
    end
    % x from G + u + omega (1 - 2x) = 0, u = ln(x/(1-x)): bisection + Newton
    a = abs(om) + 1; lo = min(max(-g - a, -700), 700); hi = min(max(-g + a, -700), 700);
    for k = 1:34
        mid = 0.5 * (lo + hi);
        pos = (g + mid + om .* (1 - 2 ./ (1 + exp(-mid)))) > 0;
        hi(pos) = mid(pos); lo(~pos) = mid(~pos);
    end
    u = 0.5 * (lo + hi);
    for k = 1:3
        xx = 1 ./ (1 + exp(-u));
        fu = g + u + om .* (1 - 2 * xx);
        u = min(max(u - fu ./ (1 - 2 * om .* xx .* (1 - xx)), -700), 700);
    end
    xx = 1 ./ (1 + exp(-u));
    qq = max(xx .* (1 - xx), 1e-300);
    ent = -(max(u, 0) + log1p(exp(-abs(u)))) + xx .* u;
    Pxx = 1 ./ qq - 2 * om; Pxd = gd + omd .* (1 - 2 * xx); Pxt = gt;
    q.p  = xx .* g + ent + om .* qq;
    q.d  = xx .* gd + omd .* qq;
    q.t  = xx .* gt;
    q.dd = xx .* gdd + omdd .* qq - Pxd.^2 ./ Pxx;
    q.tt = xx .* gtt - Pxt.^2 ./ Pxx;
    q.dt = xx .* gdt - Pxd .* Pxt ./ Pxx;
end
