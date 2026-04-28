function out = SF_getprop(PT, material, props)
% Version 1.1.0 ; Journaux et al. 2025
% Baptiste Journaux - 2026
% Calculate thermodynamic quantities for water or ice polymorphs
% (Ih, II, III, V, VI, VII and X) and for aqueous NaCl solution.
% Needs SeaFreeze_Gibbs.mat (water/ices) and SeaFreeze_Gibbs_VII_NaCl5GPa.mat
% (NaClaq) on the MATLAB path.
%
% Reference publication: Journaux et al., (2020)
%
% This is the primary entry point as of v1.1.0. The legacy name `SeaFreeze`
% is kept as a deprecated alias that forwards here with a one-time warning.
%
% Usage:
%    out = SF_getprop(PT, material)           % all supported properties
%    out = SF_getprop(PT, material, props)    % only the listed properties
%
% where:
%   out is a structure containing (SI units):
%       - G (J/kg)           Gibbs energy
%       - S (J/K/kg)         Entropy
%       - U (J/kg)           Internal energy
%       - H (J/kg)           Enthalpy
%       - A (J/kg)           Helmholtz free energy
%       - rho (kg/m^3)       Density
%       - Cp (J/kg/K)        Isobaric heat capacity
%       - Cv (J/kg/K)        Isochoric heat capacity
%       - Kt (MPa)           Isothermal bulk modulus
%       - Kp (-)             Pressure derivative of Kt
%       - Ks (MPa)           Isoentropic bulk modulus
%       - alpha (1/K)        Thermal expansivity
%       - vel (m/s)          Bulk sound speed
%     For solid phases additionally (computed only if requested):
%       - shear (MPa), Vp (m/s), Vs (m/s)
%     For NaClaq additionally:
%       - mus, muw, Va, Cpa, Vm, Cpm, phi, Vex, aw, f, m
%
%   PT is a cell {P,T} (gridded) or an N-by-2 matrix (scatter) for pure
%       phases. For NaClaq, include molality: {P,T,m} or N-by-3 [P T m].
%       P in MPa, T in K, m in mol/kg. NaN outside parametrization bounds.
%
%   props is an optional cell array of property names (or single string).
%       Omit or pass [] to compute all supported properties.
%       Example: out = SF_getprop({P,T}, 'VI', {'rho','Cp'});
%                out = SF_getprop([900 255], 'VI', 'G');
%
% Material options:
%   Ih, II, III, V, VI                  Journaux et al. 2020 / Feistel & Wagner 2006
%   VII_X_French                        French and Redmer 2015
%   water1                              Bollengier et al. 2019 (<=500 K, <=2300 MPa)
%   water2                              Brown 2018 (up to 100 GPa)
%   water_IAPWS95                       IAPWS95, Wagner & Pruss 2002
%   NaClaq                              aqueous NaCl solution (3D LBF)
%
% The ice Gibbs parametrizations are optimized for use with 'water1'
% (Bollengier et al. 2019) — use that for melting curves in 200-355 K
% and up to 2300 MPa.
%
% Examples:
%   - All properties for ice VI at 900 MPa, 255 K:
%         out = SF_getprop([900 255], 'VI')
%   - Only rho and Cp for an ice V grid:
%         out = SF_getprop({400:2:500, 240:0.5:250}, 'V', {'rho','Cp'})
%   - Only G for NaClaq:
%         out = SF_getprop({0.1:10:500, 273:400, 0.1:0.5:5}, 'NaClaq', 'G')
%
% References:
%   Bollengier, Brown & Shaw (2019) J. Chem. Phys. 151
%   Brown (2018) Fluid Phase Equilibria 463, 18-31
%   Feistel & Wagner (2006) J. Phys. Chem. Ref. Data 35, 1021-1047
%   Journaux et al. (2020) JGR: Planets 125, e2019JE006176
%   Wagner & Pruss (2002) J. Phys. Chem. Ref. Data 31, 387-535
%   French & Redmer (2015) Phys. Rev. B 91, 014308

% ---- Resolve props argument and determine what to ask the evaluator for ----
all_base_props = {'G','S','U','H','A','rho','Cp','Cv','Kt','Kp','Ks','alpha','vel','Js','gamma_Gruneisen'};
shear_props    = {'shear','Vp','Vs'};
mixing_props   = {'mus','muw','f','m','xs','xw','Va','Cpa','Vm','Vw','Cpm','phi','Vex','aw'};
solid_phases   = {'Ih','II','III','V','VI','VII_X_French'};
known_materials = [solid_phases, {'water1','water2','water_IAPWS95','NaClaq'}];

% --- Validate material -----------------------------------------------------
if ~(ischar(material) || (isstring(material) && isscalar(material)))
    error('SeaFreeze:badInput', '''material'' must be a string or character vector.');
end
material = char(material);
if ~ismember(material, known_materials)
    error('SeaFreeze:unknownMaterial', ...
          'Unknown material ''%s''. Valid options: %s', ...
          material, strjoin(known_materials, ', '));
end

% --- Validate PT shape & contents ------------------------------------------
sf_validate_PT(PT, material);

% --- Resolve and validate props --------------------------------------------
if nargin < 3 || isempty(props)
    user_props = {};   % empty => request all (and evaluator will return all)
else
    if ischar(props) || isstring(props)
        user_props = cellstr(props);
    elseif iscell(props)
        user_props = props;
    else
        error('SeaFreeze:badInput', ...
              '''props'' must be a string, char vector, or cell array of names.');
    end
end

% Deprecation alias: 'F' (Helmholtz) was an alias for 'A' in early v1.1.x.
% Translate it transparently and emit a one-shot warning per session.
user_wants_F = ~isempty(user_props) && any(strcmp(user_props, 'F'));
user_wants_A = ~isempty(user_props) && any(strcmp(user_props, 'A'));
if user_wants_F
    persistent warned_F %#ok<TLEV>
    if isempty(warned_F)
        warning('SeaFreeze:deprecatedProperty', ...
            ['Property name ''F'' is deprecated; use ''A'' (Helmholtz energy) instead. ' ...
             'Suppress with: warning(''off'',''SeaFreeze:deprecatedProperty'')']);
        warned_F = true;
    end
    user_props = union(setdiff(user_props, {'F'}, 'stable'), {'A'}, 'stable');
end

% Build the set of valid property names for this material.
if ismember(material, solid_phases)
    valid_props = [all_base_props, shear_props];
elseif strcmp(material, 'NaClaq')
    valid_props = [all_base_props, mixing_props];
else  % liquid waters: no shear, no mixing
    valid_props = all_base_props;
end

if ~isempty(user_props)
    if ~iscellstr(user_props) %#ok<ISCLSTR>
        error('SeaFreeze:badInput', '''props'' entries must all be strings.');
    end
    unknown = setdiff(user_props, valid_props);
    if ~isempty(unknown)
        error('SeaFreeze:unknownProperty', ...
              ['Unknown or unsupported property name(s) for material ''%s'': %s.\n' ...
               'Valid: %s'], ...
              material, strjoin(unknown, ', '), strjoin(valid_props, ', '));
    end
end

% Determine which shear-related properties the user wants (solid phases only)
user_wants_shear = ~isempty(intersect(user_props, shear_props));

% Build the prop list passed to fnGval. shear/Vp/Vs are computed here (not
% in the evaluator) but require rho and Ks internally.
if isempty(user_props)
    eval_props = [];  % all
else
    eval_props = setdiff(user_props, shear_props, 'stable');
    if user_wants_shear
        eval_props = union(eval_props, {'rho','Ks'}, 'stable');
    end
end

% ---- NaClaq: 3D LBF spline ------------------------------------------------
if strcmp(material,'NaClaq')
    load('SeaFreeze_Gibbs_VII_NaCl5GPa.mat', 'sp_NaCl_5GPa_500K');
    sp = sp_NaCl_5GPa_500K;
    if ~isfield(sp,'MW'), sp.MW = 58.44e-3; end
    if ~isfield(sp,'nu'), sp.nu = 2;        end
    if isempty(eval_props)
        out = fnGval(sp, PT);
    else
        out = fnGval(sp, PT, eval_props);
    end
    if user_wants_F
        out.F = out.A;
        if ~user_wants_A, out = rmfield(out, 'A'); end
    end
    return
end

% ---- Load spline for requested water/ice phase -----------------------------
shear_mod = [];
switch material
    case 'Ih'
        load('SeaFreeze_Gibbs.mat', 'G_iceIh');            sp = G_iceIh;
        shear_mod = [3.1  -0.00462 0        -0.00657 1000 273.15];
    case 'II'
        load('SeaFreeze_Gibbs.mat', 'G_iceII');            sp = G_iceII;
        shear_mod = [4.1   0.0175  0        -0.014   1100 273];
    case 'III'
        load('SeaFreeze_Gibbs.mat', 'G_iceIII');           sp = G_iceIII;
        shear_mod = [2.57  0.0175  0        -0.014   1100 273];
    case 'V'
        load('SeaFreeze_Gibbs.mat', 'G_iceV');             sp = G_iceV;
        shear_mod = [2.57  0.0175  0        -0.014   1100 273];
    case 'VI'
        load('SeaFreeze_Gibbs.mat', 'G_iceVI');            sp = G_iceVI;
        shear_mod = [2.57  0.0175  0        -0.014   1100 273];
    case 'VII_X_French'
        load('SeaFreeze_Gibbs.mat', 'G_iceVII_X_French');  sp = G_iceVII_X_French;
        shear_mod = [10    0.0033  0.000048 -0.014   1300 273];
    case 'water1'
        load('SeaFreeze_Gibbs.mat', 'G_H2O_2GPa_500K');    sp = G_H2O_2GPa_500K;
    case 'water2'
        load('SeaFreeze_Gibbs.mat', 'G_H2O_100GPa_10000K');sp = G_H2O_100GPa_10000K;
    case 'water_IAPWS95'
        load('SeaFreeze_Gibbs.mat', 'G_H2O_IAPWS');        sp = G_H2O_IAPWS;
    otherwise
        % unreachable: material was validated against known_materials above
        error('SeaFreeze:unknownMaterial', 'Unknown material ''%s''', material);
end

% fnGval uses sp_val internally — no Curve Fitting Toolbox required.
if isempty(eval_props)
    out = fnGval(sp, PT);
else
    out = fnGval(sp, PT, eval_props);
end

% ---- Shear modulus / Vp / Vs for solid phases ------------------------------
if ~isempty(shear_mod) && (isempty(user_props) || user_wants_shear)
    if iscell(PT)
        [~,Tm] = ndgrid(PT{1}, PT{2});
    else
        Tm = PT(:,2);
    end
    shear = shear_mod(1) + shear_mod(2)*(out.rho - shear_mod(5)) + ...
            shear_mod(3)*(out.rho - shear_mod(5)).^2 + ...
            shear_mod(4)*(Tm - shear_mod(6));
    want = @(name) isempty(user_props) || any(strcmp(user_props, name));
    if want('Vp'),    out.Vp    = 1e3*sqrt((out.Ks/1e3 + 4/3*shear)./out.rho/1e-3); end
    if want('Vs'),    out.Vs    = 1e3*sqrt(shear./out.rho/1e-3);                     end
    if want('shear'), out.shear = shear*1e3;  end  % -> MPa

    % Strip rho/Ks if they were only added internally to support shear calcs.
    if ~isempty(user_props)
        if ~any(strcmp(user_props,'rho')) && isfield(out,'rho'), out = rmfield(out,'rho'); end
        if ~any(strcmp(user_props,'Ks'))  && isfield(out,'Ks'),  out = rmfield(out,'Ks');  end
    end
end

% Reconstitute deprecated 'F' alias if the user requested it.
if user_wants_F
    out.F = out.A;
    if ~user_wants_A, out = rmfield(out, 'A'); end
end
