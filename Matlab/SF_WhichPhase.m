function out = SF_WhichPhase(PT, varargin)
% Version 1.1.0
% Baptiste Journaux - 2026
% Determine which supported phase is stable at given conditions.
%
% Usage:
%   out = SF_WhichPhase({P,T})                       % default: pure water
%   out = SF_WhichPhase([P T])                       % scatter form
%   out = SF_WhichPhase({P,T,m}, 'solute','NaCl')    % NaCl(aq) liquid
%   out = SF_WhichPhase([P T m], 'solute','NaCl')
%
% Returns a numerical phase code:
%   0 = liquid  (pure water1 by default, NaCl(aq) when 'solute','NaCl')
%   1 = ice Ih
%   2 = ice II
%   3 = ice III
%   5 = ice V
%   6 = ice VI
%   NaN where (P,T[,m]) is outside the parametrization range of every phase.
%
% In NaCl mode the comparison is done on the chemical potential of water:
% muw_solution (J/mol of H2O) vs G_ice * MW_H2O (J/mol). This accounts for
% freezing-point depression by salt.
%
% Example:
%   % Pure water, every 10 MPa from 0 to 2300 MPa, every 1 K from 200-355 K
%   out = SF_WhichPhase({0:10:2300, 200:355});
%
%   % NaCl(aq) solution at 2 mol/kg, sweeping P,T
%   out = SF_WhichPhase({0:10:1000, 240:1:300, 2}, 'solute','NaCl');

p = inputParser;
addParameter(p, 'solute', 'none', @(s) ischar(s) || (isstring(s) && isscalar(s)));
parse(p, varargin{:});
solute = char(p.Results.solute);
if ~ismember(lower(solute), {'none','nacl','naclaq'})
    error('SF_WhichPhase:badInput', ...
          '''solute'' must be ''none'' or ''NaCl'' (got ''%s'').', solute);
end
is_nacl = strcmpi(solute, 'NaCl') || strcmpi(solute, 'NaClaq');

% Validate PT against the implied material.
if is_nacl, mat_check = 'NaClaq'; else, mat_check = 'water1'; end
sf_validate_PT(PT, mat_check);

load('SeaFreeze_Gibbs.mat', 'G_iceIh','G_iceII','G_iceIII','G_iceV','G_iceVI', ...
     'G_H2O_2GPa_500K');

MW_H2O = 0.018015268;  % kg/mol

% --- prepare (P,T) coordinates for the ice splines and evaluate the liquid -
if is_nacl
    if iscell(PT)
        if length(PT) ~= 3
            error('SF_WhichPhase:badInput', ...
                  'NaClaq mode requires PT = {P,T,m}');
        end
        PT_pt = {PT{1}(:)', PT{2}(:)'};
        is_grid = true;
        nm = numel(PT{3});
    else
        if size(PT,2) ~= 3
            error('SF_WhichPhase:badInput', ...
                  'NaClaq mode requires PT = [P T m] (N-by-3)');
        end
        PT_pt = PT(:,1:2);
        is_grid = false;
    end
    sol  = SF_getprop(PT, 'NaClaq', 'muw');
    Gliq = sol.muw;       % J/mol of H2O in solution
    ice_scale = MW_H2O;   % convert ice G (J/kg) to J/mol of H2O
else
    if iscell(PT)
        PT_pt = {PT{1}(:)', PT{2}(:)'};
        is_grid = true;
    else
        PT_pt = PT(:,1:2);
        is_grid = false;
    end
    Gliq = fnval(G_H2O_2GPa_500K, PT_pt');   % J/kg
    ice_scale = 1;
end

% --- evaluate ice phases on the (P,T) grid / scatter list -------------------
Gih  = fnval(G_iceIh,  PT_pt') * ice_scale;
Gii  = fnval(G_iceII,  PT_pt') * ice_scale;
Giii = fnval(G_iceIII, PT_pt') * ice_scale;
Gv   = fnval(G_iceV,   PT_pt') * ice_scale;
Gvi  = fnval(G_iceVI,  PT_pt') * ice_scale;

% --- broadcast ice values across the molality axis if needed ---------------
if is_nacl && is_grid
    expand = @(X) repmat(X, [1 1 nm]);
    Gih  = expand(Gih);
    Gii  = expand(Gii);
    Giii = expand(Giii);
    Gv   = expand(Gv);
    Gvi  = expand(Gvi);
end

% --- pick the phase with minimum Gibbs energy / chemical potential ---------
sz = size(Gliq);
G  = [Gliq(:) Gih(:) Gii(:) Giii(:) Gv(:) Gvi(:)];
G(G == 0) = NaN;          % preserve original out-of-range sentinel handling
phase_codes = [0 1 2 3 5 6];

[~, idx] = min(G, [], 2, 'omitnan');
out = phase_codes(idx);
out(all(isnan(G), 2)) = NaN;
out = reshape(out, sz);
end
