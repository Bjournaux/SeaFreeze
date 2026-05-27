function P_out = SF_rho2P(rho_target, T, material, varargin)
% SF_rho2P  Invert the SeaFreeze EOS: find P such that rho(P,T) == rho_target.
% Version 1.0
% Baptiste Journaux - 2026
%
% Usage:
%   P = SF_rho2P(rho, T, material)
%   P = SF_rho2P(rho, T, material, m)              % NaClaq: molality (mol/kg)
%   P = SF_rho2P(rho, T, material, 'P0', Pguess)   % initial P guess (MPa)
%   P = SF_rho2P(rho, T, material, 'tol', 0.001)   % convergence tol (MPa, default 0.01)
%   P = SF_rho2P(rho, T, material, m, 'tol', 0.001)
%
% Inputs:
%   rho_target  – target density (kg/m³), scalar or array
%   T           – temperature (K), scalar or array (broadcast against rho_target)
%   material    – any SF_getprop material string
%                 ('water1', 'water2', 'water_IAPWS95',
%                  'Ih','II','III','V','VI','VII_X_French',
%                  'NaClaq', 'NaClaq_LP', 'NaClaq_HP', 'NaClaq_5GPa_2024')
%   m           – (NaClaq only) molality in mol/kg, scalar or array matching
%                 rho_target (positional 4th argument before any name-value pairs)
%
% Name-value options:
%   'P0'     – scalar or array initial pressure guess (MPa).
%              If omitted, an estimate from one Newton step at P_lo is used.
%   'tol'    – convergence tolerance in MPa (default: 0.01 MPa).
%   'maxiter'– max Newton-Raphson iterations before bisection fallback (default: 30).
%
% Output:
%   P_out  – pressure (MPa), same size as rho_target.
%            NaN where no solution was found (rho_target outside the phase's
%            density range at T, or T outside the spline domain).
%
% Method:
%   Newton-Raphson using the isothermal bulk modulus Kt (MPa):
%     P_{n+1} = P_n + (rho_target - rho(P_n,T)) * Kt(P_n,T) / rho(P_n,T)
%   Bisection fallback for points that stall. Convergence is verified by a
%   final residual check; points whose |rho_final - rho_target| > 100*tol
%   are returned as NaN.
%
% Examples:
%   % Pure water near ambient
%   P = SF_rho2P(997, 298, 'water1')        % ≈ 0.1 MPa
%
%   % Ice VI scatter
%   P = SF_rho2P([1310 1350 1390], [255 260 265], 'VI')
%
%   % NaClaq at 1 mol/kg
%   P = SF_rho2P(1050, 300, 'NaClaq', 1.0)
%
% See also: SF_getprop, SF_phase_range.

% ---- Validate material -------------------------------------------------------
defs = sf_material_defs();
if ~(ischar(material) || (isstring(material) && isscalar(material)))
    error('SF_rho2P:badInput', '''material'' must be a string or char vector.');
end
material = char(material);
if ~ismember(material, defs.known_materials)
    error('SF_rho2P:unknownMaterial', ...
          'Unknown material ''%s''. Valid: %s', ...
          material, strjoin(defs.known_materials, ', '));
end

is_nacl = ismember(material, defs.nacl_materials);

% ---- Extract positional molality argument and name-value pairs ---------------
m_val    = [];
nv_start = 1;
if is_nacl && ~isempty(varargin) && isnumeric(varargin{1})
    m_val    = varargin{1};
    nv_start = 2;
end

p = inputParser;
addParameter(p, 'P0',      [],   @(x) isnumeric(x) && all(isfinite(x(:))));
addParameter(p, 'tol',     0.01, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'maxiter', 30,   @(x) isnumeric(x) && isscalar(x) && x >= 1);
parse(p, varargin{nv_start:end});
P0_user  = p.Results.P0;
tol      = p.Results.tol;
max_iter = round(p.Results.maxiter);

% ---- Validate rho_target and T -----------------------------------------------
if ~isnumeric(rho_target) || ~isreal(rho_target) || isempty(rho_target)
    error('SF_rho2P:badInput', '''rho_target'' must be a non-empty real numeric array.');
end
if ~isnumeric(T) || ~isreal(T) || isempty(T)
    error('SF_rho2P:badInput', '''T'' must be a non-empty real numeric array.');
end

sz       = size(rho_target);
rho_flat = rho_target(:);
n        = numel(rho_flat);

T_flat = T(:);
if isscalar(T_flat)
    T_flat = repmat(T_flat, n, 1);
elseif numel(T_flat) ~= n
    error('SF_rho2P:badInput', ...
          'T must be scalar or the same number of elements as rho_target (%d), got %d.', ...
          n, numel(T_flat));
end

if is_nacl
    if isempty(m_val)
        error('SF_rho2P:badInput', ...
              'NaClaq material requires a molality argument m (4th positional argument).');
    end
    m_flat = m_val(:);
    if isscalar(m_flat)
        m_flat = repmat(m_flat, n, 1);
    elseif numel(m_flat) ~= n
        error('SF_rho2P:badInput', ...
              'm must be scalar or match numel(rho_target) (%d), got %d.', ...
              n, numel(m_flat));
    end
    if any(m_flat < 0)
        error('SF_rho2P:badInput', 'Molality must be non-negative.');
    end
else
    m_flat = [];
end

% ---- Domain bounds -----------------------------------------------------------
rng  = SF_phase_range(material);
P_lo = rng.P(1);
P_hi = rng.P(2);

% P_newton_lo: domain floor for Newton clipping and bisection.
% A tiny positive epsilon is used when P_lo == 0 to avoid the very few EOS
% (e.g. water2) that diverge at P = 0 exactly.
P_newton_lo = P_lo;
if P_lo == 0
    P_newton_lo = 1e-3;   % 0.001 MPa — well within any supported domain
end

% P_search_lo: interior starting point for the initial Newton-step guess only.
% For wide domains (water2: P_hi = 100 000 MPa) this keeps the seed away from
% the divergent low-P tail; for narrow ones it equals P_newton_lo.
P_search_lo = max(P_newton_lo, (P_hi - P_lo) * 0.001);

% ---- Initial guess -----------------------------------------------------------
if ~isempty(P0_user)
    P0    = P0_user(:);
    if isscalar(P0), P0 = repmat(P0, n, 1); end
    P_n   = max(P_newton_lo, min(P_hi, P0));
else
    % One linearised Newton step from a stable interior point gives a good guess.
    P_start = repmat(P_search_lo, n, 1);
    [rho0, Kt0] = eval_rho_Kt(P_start, T_flat, material, m_flat);
    dP0 = (rho_flat - rho0) .* Kt0 ./ max(abs(rho0), 1);
    P_n = P_search_lo + dP0;
    P_n = max(P_search_lo, min(P_hi, P_n));
    bad = ~isfinite(P_n);
    P_n(bad) = (P_search_lo + P_hi) / 2;
end

% ---- Newton-Raphson ----------------------------------------------------------
converged = false(n, 1);
for iter = 1:max_iter
    [rho_n, Kt_n] = eval_rho_Kt(P_n, T_flat, material, m_flat);
    residual = rho_flat - rho_n;
    dP       = residual .* Kt_n ./ max(abs(rho_n), 1);
    dP(~isfinite(dP)) = 0;
    P_new    = P_n + dP;
    P_new    = max(P_newton_lo, min(P_hi, P_new));
    converged = converged | (abs(dP) < tol);
    P_n      = P_new;
    if all(converged), break; end
end

% ---- Bisection fallback for non-converged points -----------------------------
need_bisect = ~converged & isfinite(rho_flat);
if any(need_bisect)
    idx    = find(need_bisect);
    P_n(idx) = bisect_solve(rho_flat(idx), T_flat(idx), material, ...
                             m_flat_at(m_flat, idx), ...
                             P_newton_lo, P_hi, tol, 60);
end

% ---- Final residual check: mark unresolved points NaN -----------------------
% Only evaluate rho at points where P_n is finite and positive.
valid_P   = isfinite(P_n) & P_n >= 0;
bad_final = true(n, 1);                  % start pessimistic
if any(valid_P)
    idx_v          = find(valid_P);
    P_eval         = P_n(idx_v);
    rho_final      = eval_rho(P_eval, T_flat(idx_v), material, m_flat_at(m_flat, idx_v));
    close_enough   = isfinite(rho_final) & ...
                     (abs(rho_final - rho_flat(idx_v)) <= 100 * tol);
    bad_final(idx_v) = ~close_enough;
end
% Points with non-finite rho_target were never solvable.
bad_final(~isfinite(rho_flat)) = true;

P_out       = P_n;
P_out(bad_final) = NaN;
P_out       = reshape(P_out, sz);
end

% =============================================================================
% Local helpers
% =============================================================================

function [rho, Kt] = eval_rho_Kt(P_vec, T_vec, material, m_vec)
% Evaluate rho (kg/m³) and Kt (MPa) at scattered (P,T[,m]) points.
PT  = build_PT(P_vec, T_vec, material, m_vec);
out = SF_getprop(PT, material, {'rho','Kt'});
rho = out.rho(:);
Kt  = out.Kt(:);
end

function rho = eval_rho(P_vec, T_vec, material, m_vec)
% Evaluate only rho (kg/m³).
PT  = build_PT(P_vec, T_vec, material, m_vec);
out = SF_getprop(PT, material, 'rho');
rho = out.rho(:);
end

function PT = build_PT(P_vec, T_vec, material, m_vec)
% Build the N×2 or N×3 scatter matrix for SF_getprop.
defs = sf_material_defs();
if ismember(material, defs.nacl_materials)
    PT = [P_vec(:), T_vec(:), m_vec(:)];
else
    PT = [P_vec(:), T_vec(:)];
end
end

function m_out = m_flat_at(m_flat, idx)
% Safe index into m_flat (may be empty for pure phases).
if isempty(m_flat)
    m_out = [];
else
    m_out = m_flat(idx);
end
end

function P_sol = bisect_solve(rho_tgt, T_vec, material, m_vec, P_lo, P_hi, tol, max_iter)
% Bisection on f(P) = rho(P,T) - rho_tgt = 0 for a batch of scatter points.
n  = numel(rho_tgt);
a  = repmat(P_lo, n, 1);
b  = repmat(P_hi, n, 1);

fa = eval_rho(a, T_vec, material, m_vec) - rho_tgt;
fb = eval_rho(b, T_vec, material, m_vec) - rho_tgt;

no_bracket = ~isfinite(fa) | ~isfinite(fb) | (fa .* fb > 0);
if any(no_bracket)
    warning('SF_rho2P:noBracket', ...
            '%d point(s) could not be bracketed for bisection.', sum(no_bracket));
end

for iter = 1:max_iter
    c  = (a + b) / 2;
    fc = eval_rho(c, T_vec, material, m_vec) - rho_tgt;
    left          = (fa .* fc <= 0);
    b(left)       = c(left);   fb(left)  = fc(left);
    a(~left)      = c(~left);  fa(~left) = fc(~left);
    if max(abs(b - a)) < tol, break; end
end

P_sol = (a + b) / 2;
P_sol(no_bracket) = NaN;
end
