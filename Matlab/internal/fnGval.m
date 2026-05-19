function Results = fnGval(sp, input, props)
% FNGVALS  Evaluate thermodynamic properties from a Gibbs energy B-spline.
%   Selective-output version: pass a props cell array to request only the
%   properties you need, avoiding unnecessary spline derivative evaluations.
%
%   Results = fnGvals(sp, PTm)
%   Results = fnGvals(sp, PTm, props)
%
%   Computes thermodynamic properties at specified P-T or P-T-m points from
%   a tensor B-spline representation of the Gibbs energy (G).  Units are
%   mostly MKS with P in MPa; exceptions are noted per field.
%
%   INPUT
%     sp    - B-spline structure (from sp_val / spline toolbox conventions).
%             Required fields: knots, coefs, order, number, dim.
%             Optional fields:
%               Tc      - critical temperature (K); if present, the spline is
%                         parameterised in tau = log(T/Tc) rather than T.
%               MW      - solute molecular weight (g/mol); required for 3-D.
%                         If MW is a 2-element vector, MW(2) is used (solute).
%               nu      - number of ions per formula unit; required for 3-D.
%               cutoff  - reference concentration (mol/kg) for apparent
%                         property calculation; default 2e-4.
%               mask    - NaN mask array; if present, Results.mask is returned.
%               PTm_mask - cell of {P,T} or {P,T,m} grid for the mask.
%
%     PTm   - evaluation points, in one of two formats:
%               Gridded:   cell array {P, T} or {P, T, m}  (vectors)
%               Scattered: matrix [P(:), T(:)] or [P(:), T(:), m(:)]
%
%     props - (optional) cell array of property names, or a single string.
%             Omit or pass [] to request all supported properties.
%             Unknown names raise an error.
%
%   OUTPUT
%     Results - struct with the requested fields (all per kg of solution
%               unless noted):
%
%       Always available (2-D and 3-D):
%         P      Input pressures                  (MPa)
%         T      Input temperatures               (K)
%         G      Gibbs energy                     (J/kg)
%         S      Entropy                          (J/kg/K)
%         U      Internal energy                  (J/kg)
%         H      Enthalpy                         (J/kg)
%         A      Helmholtz energy                 (J/kg)
%         rho    Density                          (kg/m^3)
%         Cp     Isobaric heat capacity           (J/kg/K)
%         Cv     Isochoric heat capacity          (J/kg/K)
%         Kt     Isothermal bulk modulus          (MPa)
%         Ks     Isentropic bulk modulus          (MPa)
%         Kp     dKt/dP                           (dimensionless)
%         alpha  Thermal expansivity              (1/K)
%         vel    Sound speed                      (m/s)
%         Js     dT/dP|adiabatic                  (K/MPa)
%         gamma_Gruneisen  Gruneisen parameter    (dimensionless)
%                  = alpha * Kt / (rho * Cv)  [Kt in Pa via 1e6 factor]
%
%       Solution properties (3-D only, per mol unless noted):
%         m      Molality                         (mol/kg)
%         xw     Mole fraction of solvent         (dimensionless)
%         xs     Mole fraction of solute          (dimensionless)
%         f      kg-solution / kg-water conversion factor
%         mus    Chemical potential, solute       (J/mol)
%         muw    Chemical potential, solvent      (J/mol)
%         Vm     Partial molal volume, solute     (cm^3/mol)
%         Vw     Partial molal volume, solvent    (cm^3/mol)
%         Va     Apparent molar volume            (cm^3/mol)
%         Vex    Excess volume                    (cm^3/mol)
%         Cpa    Apparent molar heat capacity     (J/K/mol)
%         Cpm    Partial molal heat capacity      (J/K/mol)
%         phi    Osmotic coefficient              (dimensionless)
%         aw     Solvent activity                 (dimensionless)
%
%       Optional:
%         mask   NaN mask interpolated to evaluation points.
%
%   NOTES ON PARTIAL MOLAL PROPERTIES
%     Apparent quantities (Va, Cpa) are numerically well-behaved because
%     they require only a difference between the pure-solvent and dilute
%     solution values.  Partial molal quantities (Vm, Cpm) involve
%     concentration derivatives and become ill-conditioned as m -> 0.
%     The code approximates infinite dilution at sp.cutoff mol/kg (default
%     2e-4).  This threshold warrants periodic review against the knot spacing.
%
%   JMB 2015-2026.

    % ------------------------------------------------------------------
    % Physical constants
    % ------------------------------------------------------------------
    nw = 1000 / 18.015268;   % mol H2O per kg H2O (Tilner-Roth & Friend, NIST)
    R  = 8.3144;             % universal gas constant (J/mol/K)
    MPa2Pa = 1e6;            % unit conversion: MPa -> Pa

    % ------------------------------------------------------------------
    % Spline flags
    % ------------------------------------------------------------------
    mu_flg = length(sp.knots) > 2;   % true if 3-D (P-T-m) spline

    flgTc = isfield(sp, 'Tc');
    if flgTc
        Tc = sp.Tc;
    else
        Tc = [];
    end

    % ------------------------------------------------------------------
    % Validate 3-D spline fields
    % ------------------------------------------------------------------
    if mu_flg
        if ~isfield(sp, 'MW')
            error('fnGval:missingMW', ...
                  'A compositional spline requires sp.MW (solute molecular weight).');
        end
        M = sp.MW;
        if numel(M) == 2, M = M(2); end

        if ~isfield(sp, 'nu')
            error('fnGval:missingNu', ...
                  'A compositional spline requires sp.nu (number of disassociated ions per formula unit of solute).');
        end
        nu = sp.nu;
    else
        M  = [];
        nu = [];
    end

    % ------------------------------------------------------------------
    % Resolve requested properties -> want struct
    % ------------------------------------------------------------------
    % P and T are always returned; they live outside the gated want system
    % since they are trivially cheap and nearly always needed.
    defs = sf_material_defs();
    base_props   = defs.base_props;
    mixing_props = defs.mixing_props;

    nd_input = numel(input);
    if ~iscell(input), nd_input = size(input, 2); end

    if nd_input == 2
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
            error('fnGval:unknownProperty', ...
                  'fnGval: unsupported property name(s): %s\nValid names: %s', ...
                  strjoin(unknown, ', '), strjoin(all_props, ', '));
        end
        req = props;
    end

    want = struct();
    for ii = 1:numel(all_props),  want.(all_props{ii}) = false; end
    for ii = 1:numel(req),        want.(req{ii})        = true;  end

    % ------------------------------------------------------------------
    % Determine which mixing properties are needed
    % ------------------------------------------------------------------
    has_mixing = false;
    if nd_input == 3 && mu_flg
        for ii = 1:numel(mixing_props)
            if want.(mixing_props{ii}), has_mixing = true; break; end
        end
    end

    % ------------------------------------------------------------------
    % Dependency graph: which intermediate quantities are needed
    % ------------------------------------------------------------------
    want_chem  = has_mixing && (want.mus || want.muw || want.phi || want.aw);

    need_Cp    = want.Cp || want.Cv || want.gamma_Gruneisen || ...
                 (has_mixing && (want.Cpm || want.Cpa));
    need_rho   = want.rho || want.alpha || want.Ks || want.U || ...
                 want.A || want.H || want.gamma_Gruneisen || ...
                 (has_mixing && (want.Va || want.Vex));
    need_S     = want.S || want.U || want.A || want.H;
    need_vel   = want.vel || want.Ks;
    need_U     = want.U || want.A;
    need_alpha = want.alpha || want.gamma_Gruneisen;
    need_Kt    = want.Kt || want.Kp || want.gamma_Gruneisen;

    % Which spline derivatives are required
    need = struct();
    need.G     = want.G || need_U || want.H || want_chem;
    need.d1T   = need_S;
    need.d2T   = need_Cp || need_vel || want.Js;
    need.dPT   = need_vel || need_alpha || want.Cv || want.Js;
    need.d1P   = need_rho || need_Kt || need_vel || ...
                 (has_mixing && (want.Vm || want.Vw || want.Va || want.Vex));
    need.d2P   = need_vel || need_Kt || want.Cv;
    need.d3P   = want.Kp;
    need.d2Tm  = has_mixing && (want.Cpm || want.Cpa);
    need.dPm   = has_mixing && (want.Vm || want.Vw || want.Vex);
    need.d1Tm  = has_mixing && (want.Cpm || want.Cpa) && flgTc;
    need.dGdm  = want_chem;

    % Pack intermediate-quantity flags for calc_pure_properties
    need.Cp    = need_Cp;
    need.rho   = need_rho;
    need.S     = need_S;
    need.vel   = need_vel;
    need.U     = need_U;
    need.alpha = need_alpha;
    need.Kt    = need_Kt;

    % ------------------------------------------------------------------
    % 1. Parse input -> grids and sp_val argument
    % ------------------------------------------------------------------
    [Pm, Tm, mm, tau, tautm, tau2tm, xs, xw, ...
     nP, nT, nm, mflg, gridded, sp_arg] = ...
        parse_input(input, mu_flg, flgTc, Tc, nw, nu, sp);

    % ------------------------------------------------------------------
    % 2. Evaluate required spline derivatives (gated by need struct)
    % ------------------------------------------------------------------
    derivs = get_derivatives(sp, sp_arg, mu_flg, need);

    % ------------------------------------------------------------------
    % 3. Apply chain rule d/dT = (d/dtau)*(dtau/dT) when flgTc is set
    % ------------------------------------------------------------------
    derivs = apply_Tc_chain_rule(derivs, tautm, tau2tm, flgTc, mu_flg, need);

    % ------------------------------------------------------------------
    % 4. Compute pure thermodynamic properties (gated by want/need)
    % ------------------------------------------------------------------
    pure = calc_pure_properties(derivs, Pm, Tm, MPa2Pa, want, need);

    % ------------------------------------------------------------------
    % 5. Compute solution properties (3-D only, gated by has_mixing)
    % ------------------------------------------------------------------
    sol = struct();
    if has_mixing
        sol = calc_solution_properties(sp, derivs, Pm, Tm, tau, mm, ...
                                       pure, gridded, M, nu, nw, R, ...
                                       want, want_chem, flgTc);
    end

    % ------------------------------------------------------------------
    % 6. Assemble Results struct (only requested fields)
    % ------------------------------------------------------------------
    Results = struct();

    if want.G,               Results.G               = pure.G;               end
    if want.S,               Results.S               = pure.S;               end
    if want.U,               Results.U               = pure.U;               end
    if want.H,               Results.H               = pure.H;               end
    if want.A,               Results.A               = pure.A;               end
    if want.rho,             Results.rho             = pure.rho;             end
    if want.Cp,              Results.Cp              = pure.Cp;              end
    if want.Cv,              Results.Cv              = pure.Cv;              end
    if want.Kt,              Results.Kt              = pure.Kt;              end
    if want.Ks,              Results.Ks              = pure.Ks;              end
    if want.Kp,              Results.Kp              = pure.Kp;              end
    if want.alpha,           Results.alpha           = pure.alpha;           end
    if want.vel,             Results.vel             = pure.vel;             end
    if want.Js,              Results.Js              = pure.Js;              end
    if want.gamma_Gruneisen, Results.gamma_Gruneisen = pure.gamma_Gruneisen; end

    % P and T: return the physical grid arrays (after any tau transform)
    if want.P
        if gridded && nd_input == 2
            Results.P = input{1}(:);   % return as column vector for gridded
        elseif gridded && nd_input == 3
            Results.P = input{1}(:);
        else
            Results.P = Pm;
        end
    end
    if want.T
        if gridded && nd_input == 2
            Results.T = input{2}(:);
        elseif gridded && nd_input == 3
            Results.T = input{2}(:);
        else
            Results.T = Tm;
        end
    end

    if mu_flg && nd_input == 3
        if want.m,   Results.m   = mm;      end
        if want.xs,  Results.xs  = xs;      end
        if want.xw,  Results.xw  = xw;      end
        if want.f,   Results.f   = sol.f;   end
        if want.mus, Results.mus = sol.mus;  end
        if want.muw, Results.muw = sol.muw;  end
        if want.Vm,  Results.Vm  = sol.Vm;  end
        if want.Vw,  Results.Vw  = sol.Vw;  end
        if want.Va,  Results.Va  = sol.Va;  end
        if want.Vex, Results.Vex = sol.Vex; end
        if want.Cpa, Results.Cpa = sol.Cpa; end
        if want.Cpm, Results.Cpm = sol.Cpm; end
        if want.phi, Results.phi = sol.phi; end
        if want.aw,  Results.aw  = sol.aw;  end
    end

    % ------------------------------------------------------------------
    % 7. Strip the prepended zero-concentration padding slice (gridded 3-D)
    % ------------------------------------------------------------------
    if mflg
        Results = trim_padding(Results, mu_flg);
        mm = mm(:,:,2:end);
        Pm = Pm(:,:,2:end);
        Tm = Tm(:,:,2:end);
        nm = nm - 1;
        if want.m,  Results.m  = mm; end
        % xs and xw were computed before prepending so they are already
        % the right length — no trimming needed.
    end

    % ------------------------------------------------------------------
    % 8. Interpolate mask if provided
    % ------------------------------------------------------------------
    if isfield(sp, 'mask')
        Results.mask = apply_mask(sp, Pm, Tm, mm, gridded, nP, nT, nm);
    end

end


% ======================================================================
%  SUBFUNCTIONS
% ======================================================================

function [Pm, Tm, mm, tau, tautm, tau2tm, xs, xw, ...
          nP, nT, nm, mflg, gridded, sp_arg] = ...
         parse_input(input, mu_flg, flgTc, Tc, nw, nu, sp)
% PARSE_INPUT  Unpack the raw PTm input and build the sp_val argument.
%
%   Handles 2-D gridded, 2-D scattered, 3-D gridded, 3-D scattered.
%   sp is passed to read sp.cutoff for the apparent property reference.

    gridded = iscell(input);

    if gridded
        nVars = numel(input);
    else
        nVars = size(input, 2);
    end

    if nVars == 1
        error('fnGval:badInput', ...
              ['Input must have at least 2 columns/cells [P, T]. ' ...
               'Check that input is column-centric: [P(:), T(:), m(:)].']);
    end

    % Initialise optional outputs
    mm = []; xs = []; xw = [];
    nP = 0; nT = 0; nm = 0;
    mflg = false;
    tautm = []; tau2tm = [];

    if nVars == 2
        % ---- 2-D: P and T only ----
        if gridded
            P = input{1};  T = input{2};
            nP = numel(P); nT = numel(T);
            [Pm, Tm] = ndgrid(P, T);
            tau = transform_T(T, flgTc, Tc);
            if flgTc
                tautm  =  Tm .^ -1;
                tau2tm = -Tm .^ -2;
            end
            sp_arg = {P, tau};
        else
            Pm = input(:,1);  Tm = input(:,2);
            tau = transform_T(Tm, flgTc, Tc);
            if flgTc
                tautm  =  Tm .^ -1;
                tau2tm = -Tm .^ -2;
            end
            sp_arg = [Pm, tau];
        end

    else
        % ---- 3-D: P, T, and concentration m ----
        % Reference concentration for apparent property calculation
        if isfield(sp, 'cutoff') && sp.cutoff > 0
            m_base = sp.cutoff;
        else
            m_base = 2e-4;
        end

        if gridded
            P = input{1};  T = input{2};  m = input{3};
            nP = numel(P); nT = numel(T);

            % Mole fractions computed before any m modification
            xs = nu * m ./ (nu * m + nw);
            xw = 1 - xs;

            % Replace exact zeros with eps to avoid division by zero
            m(m == 0) = eps;

            % Prepend reference concentration for apparent property calculation
            if mu_flg && m(1) > m_base
                m = [m_base; m(:)];
                mflg = true;
            end
            nm = numel(m);

            [Pm, Tm, mm] = ndgrid(P, T, m);
            tau = transform_T(T, flgTc, Tc);
            if flgTc
                tautm  =  Tm .^ -1;
                tau2tm = -Tm .^ -2;
            end
            sp_arg = {P, tau, m};

        else
            Pm = input(:,1);  Tm = input(:,2);  mm = input(:,3);

            xs = nu * mm ./ (nu * mm + nw);
            xw = 1 - xs;

            tau = transform_T(Tm, flgTc, Tc);
            if flgTc
                tautm  =  Tm .^ -1;
                tau2tm = -Tm .^ -2;
            end
            sp_arg = [Pm, tau, mm];
        end
    end
end


function tau = transform_T(T, flgTc, Tc)
% TRANSFORM_T  Convert temperature to spline coordinate.
    if flgTc
        tau = log(T / Tc);
    else
        tau = T;
    end
end


function derivs = get_derivatives(sp, sp_arg, mu_flg, need)
% GET_DERIVATIVES  Evaluate required spline derivatives (gated by need struct).
%
%   Only calls sp_val for derivatives flagged true in `need`.
%   sp_val accepts both cell (gridded) and matrix (scattered) sp_arg
%   transparently, so no gridded/scattered branching is needed here.
%   Unused fields are set to [] so downstream code can detect them.

    function v = ev(dv)
        v = sp_val(sp, dv, sp_arg);
    end

    if ~mu_flg
        % 2-D derivatives: [dP, dT] order
        if need.G,   derivs.G   = ev([0 0]); else, derivs.G   = []; end
        if need.d1T, derivs.d1T = ev([0 1]); else, derivs.d1T = []; end
        if need.d2T, derivs.d2T = ev([0 2]); else, derivs.d2T = []; end
        if need.dPT, derivs.dPT = ev([1 1]); else, derivs.dPT = []; end
        if need.d1P, derivs.d1P = ev([1 0]); else, derivs.d1P = []; end
        if need.d2P, derivs.d2P = ev([2 0]); else, derivs.d2P = []; end
        if need.d3P, derivs.d3P = ev([3 0]); else, derivs.d3P = []; end
        derivs.d1Tm = []; derivs.d2Tm = []; derivs.dPm = []; derivs.dGdm = [];
    else
        % 3-D derivatives: [dP, dT, dm] order
        if need.G,    derivs.G    = ev([0 0 0]); else, derivs.G    = []; end
        if need.d1T,  derivs.d1T  = ev([0 1 0]); else, derivs.d1T  = []; end
        if need.d2T,  derivs.d2T  = ev([0 2 0]); else, derivs.d2T  = []; end
        if need.dPT,  derivs.dPT  = ev([1 1 0]); else, derivs.dPT  = []; end
        if need.d1P,  derivs.d1P  = ev([1 0 0]); else, derivs.d1P  = []; end
        if need.d2P,  derivs.d2P  = ev([2 0 0]); else, derivs.d2P  = []; end
        if need.d3P,  derivs.d3P  = ev([3 0 0]); else, derivs.d3P  = []; end
        if need.d1Tm, derivs.d1Tm = ev([0 1 1]); else, derivs.d1Tm = []; end
        if need.d2Tm, derivs.d2Tm = ev([0 2 1]); else, derivs.d2Tm = []; end
        if need.dPm,  derivs.dPm  = ev([1 0 1]); else, derivs.dPm  = []; end
        if need.dGdm, derivs.dGdm = ev([0 0 1]); else, derivs.dGdm = []; end
    end
end


function derivs = apply_Tc_chain_rule(derivs, tautm, tau2tm, flgTc, mu_flg, need)
% APPLY_TC_CHAIN_RULE  Convert tau-domain derivatives to T-domain.
%
%   When the spline is parameterised in tau = log(T/Tc):
%     dG/dT    = (dG/dtau) / T
%     d2G/dT2  = (d2G/dtau2) / T^2  +  (dG/dtau) * (-1/T^2)
%
%   d2T must be corrected before d1T is overwritten.
    if ~flgTc, return, end

    if need.d2T, derivs.d2T = derivs.d2T .* tautm.^2 + derivs.d1T .* tau2tm; end
    if need.d1T, derivs.d1T = derivs.d1T .* tautm; end
    if need.dPT, derivs.dPT = derivs.dPT .* tautm; end

    if mu_flg
        if need.d2Tm, derivs.d2Tm = derivs.d2Tm .* tautm.^2 + derivs.d1Tm .* tau2tm; end
        if need.d1Tm, derivs.d1Tm = derivs.d1Tm .* tautm; end
    end
end


function pure = calc_pure_properties(derivs, Pm, Tm, MPa2Pa, want, need)
% CALC_PURE_PROPERTIES  Standard thermodynamic properties from G derivatives.
%   Returns a struct; fields are only populated when want/need flags are set.

    pure = struct();

    G   = derivs.G;   d1T = derivs.d1T; d2T = derivs.d2T;
    dPT = derivs.dPT; d1P = derivs.d1P; d2P = derivs.d2P; d3P = derivs.d3P;

    if want.G,     pure.G     = G;                                           end
    if need.S,     pure.S     = -d1T;                                        end
    if need.Cp,    pure.Cp    = -d2T .* Tm;                                  end
    if want.Cv,    pure.Cv    =  pure.Cp + Tm .* dPT.^2 ./ d2P;             end
    if need.rho,   pure.rho   =  MPa2Pa .* d1P.^-1;                         end
    if need.vel,   pure.vel   =  real(sqrt(d1P.^2 ./ (dPT.^2./d2T - d2P))); end
    if want.Ks,    pure.Ks    =  pure.rho .* pure.vel.^2 / MPa2Pa;          end
    if need.alpha, pure.alpha =  dPT ./ d1P;                                 end
    if need.Kt,    pure.Kt    = -d1P ./ d2P;                                 end
    if want.Kp,    pure.Kp    =  d1P .* d2P.^(-2) .* d3P - 1;              end
    if want.Js,    pure.Js    = -dPT ./ d2T;                                 end
    if need.U
        pure.U = G - MPa2Pa .* Pm ./ pure.rho + Tm .* pure.S;
    end
    if want.A,     pure.A     =  pure.U - Tm .* pure.S;                     end
    if want.H,     pure.H     =  G + Tm .* pure.S;                          end
    if want.gamma_Gruneisen
        % Cv is not necessarily in pure yet (may not have been requested);
        % compute it locally to avoid coupling want flags.
        Cv_g = pure.Cp + Tm .* dPT.^2 ./ d2P;
        pure.gamma_Gruneisen = MPa2Pa .* pure.alpha .* pure.Kt ./ pure.rho ./ Cv_g;
    end
end


function sol = calc_solution_properties(sp, derivs, Pm, Tm, tau, mm, ...
                                         pure, gridded, M, nu, nw, R, ...
                                         want, want_chem, flgTc)
% CALC_SOLUTION_PROPERTIES  Partial molal and apparent properties (gated).
%
% Arguments:
%   sp        - spline struct (for scattered sp_val reference evaluations)
%   derivs    - struct of spline derivatives at evaluation points
%   Pm, Tm    - pressure and temperature grids/vectors
%   tau       - temperature coordinate in spline space (T or log(T/Tc))
%   mm        - concentration grid/vector (mol/kg)
%   pure      - struct of pure-fluid properties (Cp used for Cpm, Cpa)
%   gridded   - true for grid input, false for scattered
%   M         - solute molar mass (kg/mol)
%   nu        - number of dissociated ions per formula unit
%   nw        - moles of water per kg of water
%   R         - gas constant (J/mol/K)
%   want      - struct of requested-output flags
%   want_chem - true if any chemical potential output is requested
%   flgTc     - true if spline is parameterised in tau = log(T/Tc)

    sol = struct();

    G    = derivs.G;
    d1P  = derivs.d1P;
    dPm  = derivs.dPm;
    dGdm = derivs.dGdm;
    d2Tm = derivs.d2Tm;

    need_sol_rho = want.Va || want.Vex;
    need_sol_Cp  = want.Cpa || want.Cpm;

    % Conversion factor f = 1 + M*m  (kg solution / kg water)
    need_f = want.f || want_chem || want.Vm || want.Vw || ...
             want.Cpm || want.Cpa || want.Va || want.Vex;
    if need_f
        f = 1 + M * mm;
        if want.f, sol.f = f; end
    end

    % Chemical potentials (J/mol)
    if want_chem
        mus = M * G + f .* dGdm;
        muw = G / nw - (mm / nw) .* f .* dGdm;
        if want.mus, sol.mus = mus; end
        if want.muw, sol.muw = muw; end
    end

    % Partial molal volumes (cm^3/mol).
    % d1P from a MPa-knot spline has implicit units of m^3/kg per MPa,
    % which numerically equals cm^3/g.  Multiplying by 1e6 (cm^3/kg per
    % mol-fraction scaling) gives cm^3/mol — no MPa2Pa factor needed here.
    % (MPa2Pa is only needed when converting d1P to SI density via rho=1/V.)
    if want.Vm || want.Vex
        sol.Vm = M * d1P + f .* dPm;    % cm^3/mol
    end
    if want.Vw
        sol.Vw = d1P / nw - (mm / nw) .* f .* dPm;   % cm^3/mol
    end

    % Partial molal heat capacity (J/K/mol)
    if want.Cpm
        Cp = pure.Cp;
        sol.Cpm = Cp * M - f .* d2Tm .* Tm;
    end

    if gridded
        % Reference values from the prepended (:,:,1) slice
        if need_sol_rho
            V  = d1P;           % specific volume in spline units = cm^3/g
            V0 = V(:,:,1);
        end
        if need_sol_Cp
            Cp  = pure.Cp;
            Cp0 = Cp(:,:,1);
        end
        if want.Va || want.Vex
            Va = (V .* f - repmat(V0, 1, 1, size(mm,3))) ./ mm;  % cm^3/mol
            if want.Va,  sol.Va  = Va;  end
            if want.Vex, sol.Vex = Va - sol.Vm(:,:,1); end
        end
        if want.Cpa
            sol.Cpa = (Cp .* f - repmat(Cp0, 1, 1, size(mm,3))) ./ mm;
        end
        if want.phi || want.aw
            phi = nw * (G(:,:,1)/nw - muw) ./ mm ./ Tm / R / nu;
            if want.phi, sol.phi = phi; end
            if want.aw,  sol.aw  = exp(-nu * mm .* phi / nw); end
        end

    else
        % Scattered: evaluate at reference concentration per-row
        mm0    = 2e-4 * ones(size(mm));
        sp_ref = [Pm(:), tau(:), mm0(:)];

        if want.Va || want.Vex
            d1P0 = sp_val(sp, [1 0 0], sp_ref);
            V0   = d1P0;         % reference specific volume (cm^3/g = spline units)
            V    = d1P;          % solution specific volume  (cm^3/g = spline units)
            Va   = (V .* f - V0) ./ mm;   % cm^3/mol
            if want.Va, sol.Va = Va; end
        end
        if want.Vex
            dPm0 = sp_val(sp, [1 0 1], sp_ref);
            f0   = 1 + M * mm0;
            Vmo  = M * d1P0 + f0 .* dPm0;   % Vm at reference conc (cm^3/mol)
            sol.Vex = Va - Vmo;
        end
        if want.Cpa
            d2T0 = sp_val(sp, [0 2 0], sp_ref);
            if flgTc
                d1T0       = sp_val(sp, [0 1 0], sp_ref);
                tautm_ref  =  Tm .^ -1;
                tau2tm_ref = -Tm .^ -2;
                d2T0 = d2T0 .* tautm_ref.^2 + d1T0 .* tau2tm_ref;
            end
            Cp0 = -d2T0 .* Tm;
            Cp  = pure.Cp;
            sol.Cpa = (Cp .* f - Cp0) ./ mm;
        end
        if want.phi || want.aw
            Gw  = sp_val(sp, sp_ref);
            phi = nw * (Gw / nw - muw) ./ mm ./ Tm / R / nu;
            if want.phi, sol.phi = phi; end
            if want.aw,  sol.aw  = exp(-nu * mm .* phi / nw); end
        end
    end
end


function Results = trim_padding(Results, mu_flg)
% TRIM_PADDING  Remove the prepended zero-concentration slice from all 3-D fields.
    fields_3d = {'G','S','U','H','A','rho','Cp','Cv','Kt','Ks','Kp', ...
                 'alpha','vel','Js','gamma_Gruneisen'};
    if mu_flg
        fields_3d = [fields_3d, {'Va','Cpa','mus','muw','f','Vm','Vw','Cpm','phi','Vex','aw'}];
    end
    for i = 1:numel(fields_3d)
        fn = fields_3d{i};
        if isfield(Results, fn) && ~isempty(Results.(fn))
            Results.(fn) = Results.(fn)(:,:,2:end);
        end
    end
end


function mask_out = apply_mask(sp, Pm, Tm, mm, gridded, nP, nT, nm)
% APPLY_MASK  Interpolate sp.mask to the evaluation points.
    mask = sp.mask;
    nd   = ndims(mask);

    if nd == 2
        pg = sp.PTm_mask{1};
        tg = sp.PTm_mask{2};
        mask = mask';
        if gridded
            mask_out = reshape(interp2(pg, tg, mask, Pm(:), Tm(:)), nP, nT);
        else
            mask_out = interp2(pg, tg, mask, Pm(:), Tm(:));
        end

    elseif nd == 3
        pg = sp.PTm_mask{1};
        tg = sp.PTm_mask{2};
        mg = sp.PTm_mask{3};
        mask = permute(mask, [2 1 3]);
        if gridded
            mask_out = reshape(interp3(pg, tg, mg, mask, Pm(:), Tm(:), mm(:)), nP, nT, nm);
        else
            mask_out = interp3(pg, tg, mg, mask, Pm(:), Tm(:), mm(:));
        end
    else
        error('fnGval:badMaskDims', 'sp.mask must be 2-D or 3-D.');
    end
end
