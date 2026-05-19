function Results = fnGval_1p0(sp, input, props)
% Frozen v1.0 baseline of fnGval, kept for regression-comparison tests
% against the upgraded fnGval. Do not modify.
%
% function to return thermodynamic properties for G splines in either
% (P and T) or (P, T, m). ALL MKS with P in MPa.
%
%   Usage:
%       Results = fnGval(sp, input)          % all supported properties
%       Results = fnGval(sp, input, props)   % only the listed properties
%
%   props : cell array of property names (or single char/string). Pass []
%           or omit to request all properties.
%     2D (PT) splines : G, S, U, H, A, rho, Cp, Cv, Kt, Kp, Ks, alpha, vel
%     3D (PTm) adds   : mus, muw, f, m, Va, Cpa, Vm, Cpm, gamma, phi,
%                       G_ss, Vex, Gex, aw
%
%   input : cell {P,T} or {P,T,m} for gridded output, or a npts-by-(2 or 3)
%           matrix [P T (m)] for scattered output. m in molality, P in MPa,
%           T in K.
%
%   Compositional (3D) splines must have fields MW (kg/mol), nu (# of ions
%   in solution), and Go (scalar or univariate LBF of the solute standard
%   state vs T).
%
% JMB 2015-2019

nw = 1000/18.01528;  % moles of water per kg of water
R  = 8.3144;         % J/mol/K

mu_flg = 1;
if length(sp.knots) == 2
    mu_flg = 0;
end

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
base_props   = {'G','S','U','H','A','rho','Cp','Cv','Kt','Kp','Ks','alpha','vel'};
mixing_props = {'mus','muw','f','m','Va','Cpa','Vm','Cpm','gamma','phi','G_ss','Vex','Gex','aw'};
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

% Intermediate quantities needed to evaluate requested outputs
need_Cp  = want.Cp  || want.Cv  || (has_mixing && (want.Cpm || want.Cpa));
need_rho = want.rho || want.alpha || want.Ks || want.U || want.A || ...
           (has_mixing && (want.Va || want.Vex));
need_S   = want.S || want.U || want.A || want.H;
need_vel = want.vel || want.Ks;
need_U   = want.U   || want.A;

% Required spline derivatives
need_G   = want.G || need_U || want.H || ...
           (has_mixing && (want.mus || want.muw || want.G_ss || want.gamma || ...
                            want.phi || want.aw  || want.Gex));
need_d1T = need_S;
need_d2T = need_Cp || need_vel;
need_d1P = need_rho || want.Kt || want.Kp || need_vel || ...
           (has_mixing && want.Vm);
need_dPT = need_vel || want.alpha || want.Cv;
need_d2P = need_vel || want.Kt || want.Kp || want.Cv;
need_d3P = want.Kp;
need_d3Tm = has_mixing && want.Cpm;
need_d2Pm = has_mixing && (want.Vm || want.Vex);
need_dGdm = has_mixing && (want.mus || want.muw || want.G_ss || want.gamma || ...
                            want.phi || want.aw || want.Gex);

%----- 2D (P,T) -------------------------------------------------------------
if nd == 2
    if iscell(input)  % gridded
        P = input{1}; T = input{2};
        [Pm,Tm] = ndgrid(P,T);
        x = {P,T};
        if need_G,   G   = fnval(sp, x);                end
        if need_d1T, d1T = fnval(fnder(sp,[0 1]), x);   end
        if need_d2T, d2T = fnval(fnder(sp,[0 2]), x);   end
        if need_d1P, d1P = fnval(fnder(sp,[1 0]), x);   end
        if need_dPT, dPT = fnval(fnder(sp,[1 1]), x);   end
        if need_d2P, d2P = fnval(fnder(sp,[2 0]), x);   end
        if need_d3P, d3P = fnval(fnder(sp,[3 0]), x);   end
    else              % scatter
        Pm = input(:,1); Tm = input(:,2);
        xt = input';
        if need_G,   G   = fnval(sp, xt)';                 end
        if need_d1T, d1T = fnval(fnder(sp,[0 1]), xt)';    end
        if need_d2T, d2T = fnval(fnder(sp,[0 2]), xt)';    end
        if need_d1P, d1P = fnval(fnder(sp,[1 0]), xt)';    end
        if need_dPT, dPT = fnval(fnder(sp,[1 1]), xt)';    end
        if need_d2P, d2P = fnval(fnder(sp,[2 0]), xt)';    end
        if need_d3P, d3P = fnval(fnder(sp,[3 0]), xt)';    end
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
    if ~isfield(sp,'Go')
        error('a compositional LBF needs the standard state of the solute at the reference P in the form of a univarient LBF')
    end
    Goin = sp.Go;

    if iscell(input)  % gridded
        P = input{1}; T = input{2}; m = input{3};
        if isstruct(Goin)
            Go = fnval(Goin, T);
        else
            Go = 1;
        end
        m(m==0) = eps;  % avoid divide-by-zero
        % Baseline molality for apparent-quantity calculation. Use sp.cutoff
        % if present (matches the Python convention); otherwise fall back to
        % machine eps as the original code did.
        if isfield(sp,'cutoff') && sp.cutoff > 0
            m_base = sp.cutoff;
        else
            m_base = eps;
        end
        mflg = 0;
        if mu_flg && (m(1) > m_base)
            m = [m_base; m(:)];
            mflg = 1;
        end
        [Pm,Tm,mm] = ndgrid(P,T,m);
        x = {P,T,m};
        if need_G,    G    = fnval(sp, x);                 end
        if need_d1T,  d1T  = fnval(fnder(sp,[0 1 0]), x);  end
        if need_d2T,  d2T  = fnval(fnder(sp,[0 2 0]), x);  end
        if need_d3Tm, d3Tm = fnval(fnder(sp,[0 2 1]), x);  end
        if need_d1P,  d1P  = fnval(fnder(sp,[1 0 0]), x);  end
        if need_d2Pm, d2Pm = fnval(fnder(sp,[1 0 1]), x);  end
        if need_dPT,  dPT  = fnval(fnder(sp,[1 1 0]), x);  end
        if need_d2P,  d2P  = fnval(fnder(sp,[2 0 0]), x);  end
        if need_d3P,  d3P  = fnval(fnder(sp,[3 0 0]), x);  end
        if need_dGdm, dGdm = fnval(fnder(sp,[0 0 1]), x);  end
    else              % scatter
        Pm = input(:,1); Tm = input(:,2); mm = input(:,3);
        if isstruct(Goin)
            Go = fnval(Goin, Tm);
        else
            Go = 1;
        end
        m = mm;
        if mu_flg
            mm0 = zeros(size(mm)) + eps;
            mflg = 1;
        else
            mflg = 0;
        end
        xt = input';
        if need_G,    G    = fnval(sp, xt)';                   end
        if need_d1T,  d1T  = fnval(fnder(sp,[0 1 0]), xt)';    end
        if need_d2T,  d2T  = fnval(fnder(sp,[0 2 0]), xt)';    end
        if need_d3Tm, d3Tm = fnval(fnder(sp,[0 2 1]), xt)';    end
        if need_d1P,  d1P  = fnval(fnder(sp,[1 0 0]), xt)';    end
        if need_d2Pm, d2Pm = fnval(fnder(sp,[1 0 1]), xt)';    end
        if need_dPT,  dPT  = fnval(fnder(sp,[1 1 0]), xt)';    end
        if need_d2P,  d2P  = fnval(fnder(sp,[2 0 0]), xt)';    end
        if need_d3P,  d3P  = fnval(fnder(sp,[3 0 0]), xt)';    end
        if need_dGdm, dGdm = fnval(fnder(sp,[0 0 1]), xt)';    end
        % Zero-concentration baselines for apparent volume / Cp (scatter)
        if mu_flg && has_mixing && (want.Va || want.Vex || want.Cpa)
            base_in = [Pm(:) Tm(:) mm0(:)]';
            if want.Cpa
                d2T0 = fnval(fnder(sp,[0 2 0]), base_in)';
                Cp0  = -d2T0.*Tm(:);
            end
            if want.Va || want.Vex
                d1P0 = fnval(fnder(sp,[1 0 0]), base_in)';
                V0   = 1e-6*d1P0;   % Pa -> MPa
            end
        end
    end
end

%----- Base thermodynamic quantities ----------------------------------------
if need_Cp,  Cp  = -d2T.*Tm;                                       end
if want.Cv,  Cv  =  Cp + Tm.*dPT.^2./d2P;                          end
if need_S,   S   = -d1T;                                           end
if need_vel, vel = real(sqrt(d1P.^2./(dPT.^2./d2T - d2P)));        end
if need_rho, rho = 1e6*d1P.^(-1);                                  end
if want.Ks,    Ks    = rho.*vel.^2/1e6;                            end
if want.alpha, alpha = 1e-6*dPT.*rho;                              end
if need_U,     U     = G - 1e6*Pm./rho + Tm.*S;                    end
if want.A,     A     = U - Tm.*S;                                  end
if want.H,     H     = G + Tm.*S;                                  end
if want.Kt,    Kt    = -d1P./d2P;                                  end
if want.Kp,    Kp    =  d1P.*d2P.^(-2).*d3P - 1;                   end

%----- 3D mixing block ------------------------------------------------------
if has_mixing
    want_chem = want.mus || want.muw || want.G_ss || want.gamma || ...
                want.phi || want.aw  || want.Gex;

    if iscell(input)  % gridded
        if want.f || want_chem || want.Vm || want.Cpm || want.Cpa || want.Va || want.Vex
            f = 1 + M*mm;
        end
        if want_chem
            mus   = M*G + f.*dGdm;
            G_ss  = Go + mus(:,:,1) - mus(1,:,1);
            muw   = G/nw - 1/nw*f.*mm.*dGdm;
            gamma = exp(1/R*1/nu*(mus - G_ss)./Tm - log(mm));
            phi   = -nw*(muw - G(:,:,1)/nw)./mm./Tm/R/nu;
            aw    = exp(-mm.*phi*2/nw);
            Gex   = nu*R*Tm.*mm.*(log(gamma) - (phi-1));
        end
        if want.Vm || want.Vex
            Vm = M*d1P + f.*d2Pm;
        end
        if want.Cpm
            Cpm = M*Cp - f.*d3Tm.*Tm;
        end
        if want.Cpa
            Cpa = (Cp.*f - repmat(squeeze(Cp(:,:,1)),1,1,length(m)))./mm;
        end
        if want.Va || want.Vex
            V  = rho.^-1;
            Va = 1e6*(V.*f - repmat(squeeze(V(:,:,1)),1,1,length(m)))./mm;
        end
        if want.Vex
            Vex = Va - Vm(:,:,1);
        end

        if mflg  % strip the prepended zero-concentration slice
            mm  = mm(:,:,2:end);
            if want.G,    G    = G(:,:,2:end);       end
            if want.rho,  rho  = rho(:,:,2:end);     end
            if need_Cp,   Cp   = Cp(:,:,2:end);      end
            if want.Cv,   Cv   = Cv(:,:,2:end);      end
            if need_S,    S    = S(:,:,2:end);       end
            if need_vel,  vel  = vel(:,:,2:end);     end
            if want.Ks,   Ks   = Ks(:,:,2:end);      end
            if want.alpha,alpha= alpha(:,:,2:end);   end
            if need_U,    U    = U(:,:,2:end);       end
            if want.A,    A    = A(:,:,2:end);       end
            if want.H,    H    = H(:,:,2:end);       end
            if want.Kt,   Kt   = Kt(:,:,2:end);      end
            if want.Kp,   Kp   = Kp(:,:,2:end);      end
            if want_chem
                mus   = mus(:,:,2:end);
                muw   = muw(:,:,2:end);
                G_ss  = G_ss(:,:,2:end);
                gamma = gamma(:,:,2:end);
                phi   = phi(:,:,2:end);
                aw    = aw(:,:,2:end);
                Gex   = Gex(:,:,2:end);
            end
            if want.f || want_chem || want.Vm || want.Cpm || want.Cpa || want.Va || want.Vex
                f = f(:,:,2:end);
            end
            if want.Vm || want.Vex, Vm  = Vm(:,:,2:end);  end
            if want.Cpm,            Cpm = Cpm(:,:,2:end); end
            if want.Cpa,            Cpa = Cpa(:,:,2:end); end
            if want.Va || want.Vex, Va  = Va(:,:,2:end);  end
            if want.Vex,            Vex = Vex(:,:,2:end); end
        end

    else  % scatter (mixing quantities here are flagged as broken in the original code)
        if want.f || want_chem || want.Vm || want.Cpm || want.Cpa || want.Va || want.Vex
            f = 1 + M*mm;
        end
        if want_chem
            warning('The mixing quantities (Vex, Gex, gamma, etc.) are currently broken for scatter input - needs code revision')
            mus   = M*G + f.*dGdm;
            muw   = G/nw - 1/nw*f.*mm.*dGdm;
            G_ss  = Go + R*nu*Tm.*log(mm);
            gamma = exp(1/R*1/nu*(mus - G_ss)./Tm - log(mm));
            phi   = -nw*(muw - G(:,:,1)/nw)./mm./Tm/R/nu;
            aw    = exp(-mm.*phi*2/nw);
            Gex   = nu*R*Tm.*(log(gamma) - (phi-1));
        end
        if want.Vm || want.Vex, Vm  = M*d1P + f.*d2Pm;                       end
        if want.Cpm,            Cpm = M*Cp  - f.*d3Tm.*Tm;                   end
        if want.Cpa,            Cpa = (Cp.*f - Cp0)./mm;                     end
        if want.Va || want.Vex
            V  = rho.^-1;
            Va = 1e6*(V.*f - V0)./mm;
        end
        if want.Vex,            Vex = Va - Vm(:,:,1);                        end
    end
end

%----- Build output struct with only the requested fields -------------------
Results = struct();
if want.G,     Results.G     = G;      end
if want.S,     Results.S     = S;      end
if want.U,     Results.U     = U;      end
if want.H,     Results.H     = H;      end
if want.A,     Results.A     = A;      end
if want.rho,   Results.rho   = rho;    end
if want.Cp,    Results.Cp    = Cp;     end
if want.Cv,    Results.Cv    = Cv;     end
if want.Kt,    Results.Kt    = Kt;     end
if want.Ks,    Results.Ks    = Ks;     end
if want.Kp,    Results.Kp    = Kp;     end
if want.alpha, Results.alpha = alpha;  end
if want.vel,   Results.vel   = vel;    end
if nd == 3 && mu_flg
    if want.Va,    Results.Va    = Va;    end
    if want.Cpa,   Results.Cpa   = Cpa;   end
    if want.mus,   Results.mus   = mus;   end
    if want.muw,   Results.muw   = muw;   end
    if want.f,     Results.f     = f;     end
    if want.m,     Results.m     = mm;    end
    if want.Vm,    Results.Vm    = Vm;    end
    if want.Cpm,   Results.Cpm   = Cpm;   end
    if want.gamma, Results.gamma = gamma; end
    if want.phi,   Results.phi   = phi;   end
    if want.G_ss,  Results.G_ss  = G_ss;  end
    if want.Vex,   Results.Vex   = Vex;   end
    if want.Gex,   Results.Gex   = Gex;   end
    if want.aw,    Results.aw    = aw;    end
end
