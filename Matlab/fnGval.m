function Results = fnGval(sp, input, props)
% fnGval  Thermodynamic properties from a Gibbs LBF spline.
%
% Upgraded base for SeaFreeze v1.1.x. Compared with the v1.0 baseline kept
% as fnGval_1p0.m, this version:
%   * uses sp_val instead of fnval/fnder (no Curve Fitting Toolbox dep);
%   * supports dimensionless temperature splines via sp.Tc (tau = log(T/Tc))
%     with proper chain-rule derivative adjustments;
%   * returns extra solution properties: Vw (partial molar volume of solvent),
%     Js (Joule-Thomson coefficient -dPT/d2T), gamma_Gruneisen (Gruneisen
%     parameter), xs/xw (mole fractions);
%   * interpolates an optional sp.mask field to flag points outside the
%     fit's validity domain;
%   * fixes the NaCl scatter-input mixing-quantity branch (apparent and
%     partial properties now use a true zero-concentration baseline at every
%     scatter row).
%
% Usage:
%     Results = fnGval(sp, input)            % all supported properties
%     Results = fnGval(sp, input, props)     % only the listed properties
%
% props : cell array of property names (or single char/string). Pass [] or
%         omit to request all properties.
%   2D (P,T) splines : G,S,U,H,A,rho,Cp,Cv,Kt,Kp,Ks,alpha,vel,Js,gamma_Gruneisen
%   3D (P,T,m) adds  : mus,muw,f,m,xs,xw,Va,Cpa,Vm,Vw,Cpm,phi,Vex,aw
%
% Naming note: 'gamma_Gruneisen' is the Gruneisen parameter
% alpha*Kt/(rho*Cv). The NaCl activity coefficient and excess-mixing
% quantities (gamma, G_ss, Gex) were removed in v1.1.x — they depended on
% a solute standard-state Go that the MATLAB port never computed correctly.
%
% input : cell {P,T} or {P,T,m} for gridded output, or N-by-(2 or 3) matrix
%         [P T (m)] for scattered output. P in MPa, T in K, m in mol/kg.
%
% Compositional (3D) splines must have fields MW (kg/mol) and nu.
%
% Original by JMB 2015-2026; integration / selective-prop port for
% SeaFreeze v1.1.x by Baptiste Journaux 2026.

nw = 1000/18.015268;   % moles of water per kg of water (Tillner-Roth & Friend)
R  = 8.3144;           % J/mol/K

mu_flg = length(sp.knots) > 2;

flgTc = isfield(sp,'Tc');
if flgTc, Tc = sp.Tc; end

% Determine input dimensionality
if iscell(input)
    nd = length(input);
else
    [~,nd] = size(input);
end
if nd == 1
    error('Check input PTm: fnGval finds one column of input. Input needs to be column centric: [P(:) T(:) optional m(:)]')
end

%----- Resolve requested properties -----------------------------------------
base_props   = {'G','S','U','H','A','rho','Cp','Cv','Kt','Kp','Ks','alpha','vel','Js','gamma_Gruneisen'};
mixing_props = {'mus','muw','f','m','xs','xw','Va','Cpa','Vm','Vw','Cpm','phi','Vex','aw'};
if nd == 2
    all_props = base_props;
else
    all_props = [base_props, mixing_props];
end

if nargin < 3 || isempty(props)
    req = all_props;
else
    if ischar(props) || isstring(props)
        props = cellstr(props);
    end
    unknown = setdiff(props, all_props);
    if ~isempty(unknown)
        error('fnGval: unsupported property name(s): %s', strjoin(unknown, ', '));
    end
    req = props;
end
want = struct();
for ii = 1:length(all_props), want.(all_props{ii}) = false; end
for ii = 1:length(req),       want.(req{ii})       = true;  end

% Any mixing quantity requested?
has_mixing = false;
if nd == 3 && mu_flg
    for ii = 1:length(mixing_props)
        if want.(mixing_props{ii}), has_mixing = true; break; end
    end
end

% Intermediate quantities required for the requested outputs.
need_Cp  = want.Cp || want.Cv || want.gamma_Gruneisen || ...
           (has_mixing && (want.Cpm || want.Cpa));
need_rho = want.rho || want.alpha || want.Ks || want.U || want.A || ...
           want.H || want.gamma_Gruneisen || ...
           (has_mixing && (want.Va || want.Vex));
need_S   = want.S || want.U || want.A || want.H;
need_vel = want.vel || want.Ks;
need_U   = want.U || want.A;
need_alpha = want.alpha || want.gamma_Gruneisen;
need_Kt    = want.Kt || want.Kp || want.gamma_Gruneisen;

% Required spline derivatives.
need_G   = want.G || need_U || want.H || ...
           (has_mixing && (want.mus || want.muw || want.phi || want.aw));
need_d1T = need_S;
need_d2T = need_Cp || need_vel || want.Js;
need_dPT = need_vel || need_alpha || want.Cv || want.Js;
need_d1P = need_rho || need_Kt || need_vel || ...
           (has_mixing && (want.Vm || want.Vw));
need_d2P = need_vel || need_Kt || want.Cv;
need_d3P = want.Kp;
need_d2Tm = has_mixing && want.Cpm;
need_dPm  = has_mixing && (want.Vm || want.Vw || want.Vex);
need_d1Tm = has_mixing && want.Cpm && flgTc;  % chain-rule needs d1Tm if Tc reduction is in use
need_dGdm = has_mixing && (want.mus || want.muw || want.phi || want.aw);

%----- 2D (P,T) -------------------------------------------------------------
if nd == 2
    if iscell(input)  % gridded
        P = input{1}; T = input{2};
        nP = length(P); nT = length(T);
        [Pm,Tm] = ndgrid(P,T);
        if flgTc
            tau   = log(T/Tc);
            tautm  = Tm.^-1;
            tau2tm = -Tm.^-2;
        else
            tau = T;
        end
        x_eval = {P, tau};
        if need_G,   G   = sp_val(sp, x_eval);                 end
        if need_d1T, d1T = sp_val(sp, [0 1], x_eval);          end
        if need_d2T, d2T = sp_val(sp, [0 2], x_eval);          end
        if need_d1P, d1P = sp_val(sp, [1 0], x_eval);          end
        if need_dPT, dPT = sp_val(sp, [1 1], x_eval);          end
        if need_d2P, d2P = sp_val(sp, [2 0], x_eval);          end
        if need_d3P, d3P = sp_val(sp, [3 0], x_eval);          end
    else              % scatter
        P = input(:,1); T = input(:,2);
        Pm = P; Tm = T;
        if flgTc
            tau    = log(Tm/Tc);
            tautm  = Tm.^-1;
            tau2tm = -Tm.^-2;
        else
            tau = Tm;
        end
        in_eval = [Pm tau];
        if need_G,   G   = sp_val(sp, in_eval);                end
        if need_d1T, d1T = sp_val(sp, [0 1], in_eval);         end
        if need_d2T, d2T = sp_val(sp, [0 2], in_eval);         end
        if need_d1P, d1P = sp_val(sp, [1 0], in_eval);         end
        if need_dPT, dPT = sp_val(sp, [1 1], in_eval);         end
        if need_d2P, d2P = sp_val(sp, [2 0], in_eval);         end
        if need_d3P, d3P = sp_val(sp, [3 0], in_eval);         end
    end

%----- 3D (P,T,m) -----------------------------------------------------------
elseif nd == 3
    if ~isfield(sp,'MW')
        error('a compositional LBF needs the molecular weight set in the structure')
    end
    M = sp.MW;
    if length(M) == 2, M = M(2); end
    if ~isfield(sp,'nu')
        error('a compositional LBF needs the number of ions in solution (nu) set in the structure')
    end
    nu = sp.nu;

    if iscell(input)  % gridded
        P = input{1}; T = input{2}; m = input{3};
        nP = length(P); nT = length(T); nm0 = length(m);  %#ok<NASGU>

        m(m==0) = eps;
        if isfield(sp,'cutoff') && sp.cutoff > 0
            m_base = sp.cutoff;
        else
            m_base = 2e-4;        % matches Python convention
        end
        mflg = 0;
        if mu_flg && (m(1) > m_base)
            m = [m_base; m(:)];
            mflg = 1;
        end
        [Pm,Tm,mm] = ndgrid(P,T,m);
        if flgTc
            tau    = log(T/Tc);
            tautm  = Tm.^-1;
            tau2tm = -Tm.^-2;
        else
            tau = T;
        end
        if has_mixing
            xs = nu*mm./(nu*mm + nw);
            xw = 1 - xs;
        end
        x_eval = {P, tau, m};
        if need_G,    G    = sp_val(sp, x_eval);              end
        if need_d1T,  d1T  = sp_val(sp, [0 1 0], x_eval);     end
        if need_d2T,  d2T  = sp_val(sp, [0 2 0], x_eval);     end
        if need_d2Tm, d2Tm = sp_val(sp, [0 2 1], x_eval);     end
        if need_d1Tm, d1Tm = sp_val(sp, [0 1 1], x_eval);     end
        if need_d1P,  d1P  = sp_val(sp, [1 0 0], x_eval);     end
        if need_dPm,  dPm  = sp_val(sp, [1 0 1], x_eval);     end
        if need_dPT,  dPT  = sp_val(sp, [1 1 0], x_eval);     end
        if need_d2P,  d2P  = sp_val(sp, [2 0 0], x_eval);     end
        if need_d3P,  d3P  = sp_val(sp, [3 0 0], x_eval);     end
        if need_dGdm, dGdm = sp_val(sp, [0 0 1], x_eval);     end
    else              % scatter
        Pm = input(:,1); Tm = input(:,2); mm = input(:,3);
        nP = length(Pm);  %#ok<NASGU>
        m = mm;
        if flgTc
            tau    = log(Tm/Tc);
            tautm  = Tm.^-1;
            tau2tm = -Tm.^-2;
        else
            tau = Tm;
        end
        if has_mixing
            xs = nu*mm./(nu*mm + nw);
            xw = 1 - xs;
        end
        in_eval = [Pm tau mm];

        if need_G,    G    = sp_val(sp, in_eval);              end
        if need_d1T,  d1T  = sp_val(sp, [0 1 0], in_eval);     end
        if need_d2T,  d2T  = sp_val(sp, [0 2 0], in_eval);     end
        if need_d2Tm, d2Tm = sp_val(sp, [0 2 1], in_eval);     end
        if need_d1Tm, d1Tm = sp_val(sp, [0 1 1], in_eval);     end
        if need_d1P,  d1P  = sp_val(sp, [1 0 0], in_eval);     end
        if need_dPm,  dPm  = sp_val(sp, [1 0 1], in_eval);     end
        if need_dPT,  dPT  = sp_val(sp, [1 1 0], in_eval);     end
        if need_d2P,  d2P  = sp_val(sp, [2 0 0], in_eval);     end
        if need_d3P,  d3P  = sp_val(sp, [3 0 0], in_eval);     end
        if need_dGdm, dGdm = sp_val(sp, [0 0 1], in_eval);     end

        % Zero-concentration baselines for apparent volume / Cp / chemical
        % potentials (per scatter row). Matches the Python convention.
        if mu_flg && has_mixing
            mm0 = 2e-4*ones(size(mm));
            base_in = [Pm tau mm0];
            d1P0 = sp_val(sp, [1 0 0], base_in);
            d2T0 = sp_val(sp, [0 2 0], base_in);
            d1T0 = sp_val(sp, [0 1 0], base_in);
            dPm0 = sp_val(sp, [1 0 1], base_in);
            Gw   = sp_val(sp, base_in);
            if flgTc
                d2T0 = d2T0.*tautm.^2 + d1T0.*tau2tm;
                d1T0 = d1T0.*tautm;
            end
        end
    end
end

% Apply chain-rule corrections for the dimensionless temperature variable
% tau = log(T/Tc).  d/dT = (1/T) d/dtau ; d2/dT2 = (1/T^2) d2/dtau2 - (1/T^2) d/dtau .
if flgTc
    if need_d2T, d2T = d2T.*tautm.^2 + d1T.*tau2tm; end
    if need_d1T, d1T = d1T.*tautm;                  end
    if need_dPT, dPT = dPT.*tautm;                  end
    if has_mixing
        if need_d2Tm, d2Tm = d2Tm.*tautm.^2 + d1Tm.*tau2tm; end
        if need_d1Tm, d1Tm = d1Tm.*tautm;                   end
    end
end

%----- Base thermodynamic quantities ----------------------------------------
if need_Cp,    Cp    = -d2T.*Tm;                                      end
if want.Cv,    Cv    =  Cp + Tm.*dPT.^2./d2P;                         end
if need_S,     S     = -d1T;                                          end
if need_vel,   vel   = real(sqrt(d1P.^2./(dPT.^2./d2T - d2P)));       end
if need_rho,   rho   = 1e6 ./ d1P;                                    end
if want.Ks,    Ks    = rho.*vel.^2/1e6;                               end
if need_alpha, alpha = 1e-6*dPT.*rho;                                 end
if need_U,     U     = G - 1e6*Pm./rho + Tm.*S;                       end
if want.A,     A     = U - Tm.*S;                                     end
if want.H,     H     = G + Tm.*S;                                     end
if need_Kt,    Kt    = -d1P./d2P;                                     end
if want.Kp,    Kp    =  d1P.*d2P.^(-2).*d3P - 1;                      end
if want.Js,    Js    = -dPT./d2T;                                     end
if want.gamma_Gruneisen
    Cv_g = Cp + Tm.*dPT.^2./d2P;       % always need Cv for Gruneisen
    gamma_Gruneisen = 1e6*alpha.*Kt./rho./Cv_g;
end

%----- 3D mixing block ------------------------------------------------------
if has_mixing
    want_chem = want.mus || want.muw || want.phi || want.aw;

    if iscell(input)  % gridded
        if want.f || want_chem || want.Vm || want.Vw || want.Cpm || ...
                want.Cpa || want.Va || want.Vex
            f = 1 + M*mm;
        end
        if want_chem
            mus = M*G + f.*dGdm;
            muw = G/nw - 1/nw*f.*mm.*dGdm;
            phi = -nw*(muw - G(:,:,1)/nw)./mm./Tm/R/nu;
            aw  = exp(-mm.*phi*nu/nw);
        end
        if want.Vm || want.Vex,    Vm  = M*d1P + f.*dPm;                                       end
        if want.Vw,                Vw  = d1P/nw - f.*mm.*dPm/nw;                               end
        if want.Cpm,               Cpm = M*Cp - f.*d2Tm.*Tm;                                   end
        if want.Cpa
            Cpa = (Cp.*f - repmat(squeeze(Cp(:,:,1)),1,1,length(m)))./mm;
        end
        if want.Va || want.Vex
            V  = rho.^-1;
            Va = 1e6*(V.*f - repmat(squeeze(V(:,:,1)),1,1,length(m)))./mm;
        end
        if want.Vex,               Vex = Va - Vm(:,:,1);                                       end

        if mflg  % strip the prepended zero-concentration slice
            mm = mm(:,:,2:end);
            if has_mixing && (want.xs || want.xw)
                xs = xs(:,:,2:end); xw = xw(:,:,2:end);
            end
            if want.G,             G    = G(:,:,2:end);     end
            if want.rho,           rho  = rho(:,:,2:end);   end
            if need_Cp,            Cp   = Cp(:,:,2:end);    end
            if want.Cv,            Cv   = Cv(:,:,2:end);    end
            if need_S,             S    = S(:,:,2:end);     end
            if need_vel,           vel  = vel(:,:,2:end);   end
            if want.Ks,            Ks   = Ks(:,:,2:end);    end
            if need_alpha,         alpha= alpha(:,:,2:end); end
            if need_U,             U    = U(:,:,2:end);     end
            if want.A,             A    = A(:,:,2:end);     end
            if want.H,             H    = H(:,:,2:end);     end
            if need_Kt,            Kt   = Kt(:,:,2:end);    end
            if want.Kp,            Kp   = Kp(:,:,2:end);    end
            if want.Js,            Js   = Js(:,:,2:end);    end
            if want.gamma_Gruneisen, gamma_Gruneisen = gamma_Gruneisen(:,:,2:end); end
            if want_chem
                mus = mus(:,:,2:end);
                muw = muw(:,:,2:end);
                phi = phi(:,:,2:end);
                aw  = aw(:,:,2:end);
            end
            if want.f || want_chem || want.Vm || want.Vw || want.Cpm || ...
                    want.Cpa || want.Va || want.Vex
                f = f(:,:,2:end);
            end
            if want.Vm || want.Vex, Vm  = Vm(:,:,2:end);  end
            if want.Vw,             Vw  = Vw(:,:,2:end);  end
            if want.Cpm,            Cpm = Cpm(:,:,2:end); end
            if want.Cpa,            Cpa = Cpa(:,:,2:end); end
            if want.Va || want.Vex, Va  = Va(:,:,2:end);  end
            if want.Vex,            Vex = Vex(:,:,2:end); end
        end

    else  % scatter — uses per-row zero-concentration baselines
        if want.f || want_chem || want.Vm || want.Vw || want.Cpm || ...
                want.Cpa || want.Va || want.Vex
            f = 1 + M*mm;
        end
        if want_chem
            mus = M*G + f.*dGdm;
            muw = G/nw - 1/nw*f.*mm.*dGdm;
            phi = -nw*(muw - Gw/nw)./mm./Tm/R/nu;
            aw  = exp(-mm.*phi*nu/nw);
        end
        if want.Vm || want.Vex
            Vm  = M*d1P + f.*dPm;
            Vmo = M*d1P0 + (1+M*2e-4)*dPm0;        % infinite-dilution Vm
        end
        if want.Vw,                 Vw  = d1P/nw - f.*mm.*dPm/nw;                               end
        if want.Cpm,                Cpm = M*Cp - f.*d2Tm.*Tm;                                   end
        if want.Cpa,                Cpa = (Cp.*f - (-d2T0.*Tm))./mm;                            end
        if want.Va || want.Vex
            V  = rho.^-1;
            V0 = 1e-6*d1P0;
            Va = 1e6*(V.*f - V0)./mm;
        end
        if want.Vex,                Vex = Va - Vmo;                                             end
    end
end

%----- Build output struct with only the requested fields -------------------
Results = struct();
if want.G,     Results.G     = G;     end
if want.S,     Results.S     = S;     end
if want.U,     Results.U     = U;     end
if want.H,     Results.H     = H;     end
if want.A,     Results.A     = A;     end
if want.rho,   Results.rho   = rho;   end
if want.Cp,    Results.Cp    = Cp;    end
if want.Cv,    Results.Cv    = Cv;    end
if want.Kt,    Results.Kt    = Kt;    end
if want.Ks,    Results.Ks    = Ks;    end
if want.Kp,    Results.Kp    = Kp;    end
if want.alpha, Results.alpha = alpha; end
if want.vel,   Results.vel   = vel;   end
if want.Js,    Results.Js    = Js;    end
if want.gamma_Gruneisen, Results.gamma_Gruneisen = gamma_Gruneisen; end

if nd == 3 && mu_flg
    if want.Va,    Results.Va    = Va;    end
    if want.Cpa,   Results.Cpa   = Cpa;   end
    if want.mus,   Results.mus   = mus;   end
    if want.muw,   Results.muw   = muw;   end
    if want.f,     Results.f     = f;     end
    if want.m,     Results.m     = mm;    end
    if want.xs,    Results.xs    = xs;    end
    if want.xw,    Results.xw    = xw;    end
    if want.Vm,    Results.Vm    = Vm;    end
    if want.Vw,    Results.Vw    = Vw;    end
    if want.Cpm,   Results.Cpm   = Cpm;   end
    if want.phi,   Results.phi   = phi;   end
    if want.Vex,   Results.Vex   = Vex;   end
    if want.aw,    Results.aw    = aw;    end
end

% Optional validity-domain mask interpolation -------------------------------
if isfield(sp,'mask')
    mk = sp.mask;
    d  = ndims(mk);
    if d == 2
        pg = sp.PTm_mask{1}; tg = sp.PTm_mask{2};
        if iscell(input)
            Results.mask = reshape(interp2(pg,tg,mk',Pm(:),Tm(:)), size(Pm));
        else
            Results.mask = interp2(pg,tg,mk',Pm(:),Tm(:));
        end
    elseif d == 3
        pg = sp.PTm_mask{1}; tg = sp.PTm_mask{2}; mg = sp.PTm_mask{3};
        mk = permute(mk,[2 1 3]);
        if iscell(input)
            Results.mask = reshape(interp3(pg,tg,mg,mk,Pm(:),Tm(:),mm(:)), size(Pm));
        else
            Results.mask = interp3(pg,tg,mg,mk,Pm(:),Tm(:),mm(:));
        end
    end
end
end
