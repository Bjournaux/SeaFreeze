function Results = SF_NaCl_stitch(spLP, spHP, PT, props)
% SF_NaCl_stitch  Evaluate NaCl(aq) properties from stitched LP + HP splines.
% Version 1.1.0 — Baptiste Journaux — 2026
%
% Stitches two complementary NaCl(aq) Gibbs splines that cover different P
% ranges:
%   spLP  — low-P spline,  valid for P ≤ ~1001 MPa
%   spHP  — high-P spline, valid for P ≥ ~499.9 MPa
%
% In the overlap [499.9, 1001] MPa a C∞ cosine taper blends each property
% individually:
%   w(P) = (1 + cos(π·(P−P_lo)/(P_hi−P_lo))) / 2
%   f_out = w·f_LP + (1−w)·f_HP
%
% Usage:
%   Results = SF_NaCl_stitch(spLP, spHP, {P, T, m})
%   Results = SF_NaCl_stitch(spLP, spHP, [P T m])           % N-by-3 scatter
%   Results = SF_NaCl_stitch(spLP, spHP, {P, T, m}, props)  % selective props
%
% Inputs:
%   spLP  - low-P NaCl Gibbs B-spline struct (e.g. spLP_2026)
%   spHP  - high-P NaCl Gibbs B-spline struct (e.g. spHP_2_2026_r3)
%   PT    - {P_vec, T_vec, m_vec_or_scalar} or N-by-3 [P T m] matrix
%   props - (optional) char or cell array of property names; omit for all
%
% Returns:
%   Results - struct with the same fields as fnGval; all units as in fnGval.
%             In the blend zone each property is a weighted average of the
%             two splines; outside the overlap only the valid spline is used.
%
% Supported properties (same as fnGval 3-D / NaClaq mode):
%   Base:    G S U H A rho Cp Cv Kt Kp Ks alpha vel Js gamma_Gruneisen P T
%   Mixing:  mus muw f m xs xw Va Cpa Vm Vw Cpm phi Vex aw
%
% See also: fnGval, SF_getprop

% ---- Blend-zone boundaries (MPa) ----------------------------------------
P_LO = 499.9;    % LP upper knot boundary / start of cosine taper
P_HI = 1001.0;   % HP lower knot boundary / end of cosine taper

if nargin < 4, props = []; end

% =========================================================================
% 1. Parse input
% =========================================================================
if iscell(PT)
    if numel(PT) ~= 3
        error('SF_NaCl_stitch:badInput', ...
              'Grid input must be a 3-element cell {P_vec, T_vec, m_vec}.');
    end
    P_vec   = double(PT{1}(:)');   % always row vector
    is_grid = true;
else
    if size(PT, 2) ~= 3
        error('SF_NaCl_stitch:badInput', ...
              'Scatter input must be an N-by-3 matrix [P T m].');
    end
    P_vec   = double(PT(:, 1))';   % row vector
    is_grid = false;
end
nP = numel(P_vec);

% =========================================================================
% 2. Cosine taper w(P)
% =========================================================================
w_P          = zeros(1, nP);
w_P(P_vec <= P_LO) = 1;
w_P(P_vec >= P_HI) = 0;
in_taper     = P_vec > P_LO & P_vec < P_HI;
w_P(in_taper) = 0.5 * (1 + cos(pi * (P_vec(in_taper) - P_LO) / (P_HI - P_LO)));

% Logical masks: which P values each spline is called on
need_LP = P_vec <= P_HI;    % LP is valid / evaluated here
need_HP = P_vec >= P_LO;    % HP is valid / evaluated here

% =========================================================================
% 3. Build sub-inputs and evaluate splines
% =========================================================================
if is_grid
    PT_LP = {P_vec(need_LP), PT{2}, PT{3}};
    PT_HP = {P_vec(need_HP), PT{2}, PT{3}};
else
    PT_LP = PT(need_LP, :);
    PT_HP = PT(need_HP, :);
end

has_LP = any(need_LP);
has_HP = any(need_HP);

res_LP = struct();
res_HP = struct();
if has_LP, res_LP = fnGval(spLP, PT_LP, props); end
if has_HP, res_HP = fnGval(spHP, PT_HP, props); end

% =========================================================================
% 4. Index mappings (P dimension only)
% =========================================================================
idx_LP = find(need_LP);    % rows in LP result → global P indices
idx_HP = find(need_HP);

% Blend zone, LP-only zone, HP-only zone (global P indices)
blend_global = find(need_LP & need_HP);
lp_only      = find(need_LP & ~need_HP);
hp_only      = find(~need_LP & need_HP);

% Positions within LP / HP result arrays
blend_in_lp = find(ismember(idx_LP, blend_global));
blend_in_hp = find(ismember(idx_HP, blend_global));
lp_only_in  = find(ismember(idx_LP, lp_only));
hp_only_in  = find(ismember(idx_HP, hp_only));

% Taper weights for the blend zone (column vector for broadcasting)
w_blend = w_P(blend_global)';

% =========================================================================
% 5. Assemble output struct
% =========================================================================
if has_LP
    all_fields = fieldnames(res_LP)';
else
    all_fields = fieldnames(res_HP)';
end

% P, T, m are input-coordinate echoes — handled separately below
COORD_FIELDS = {'P', 'T', 'm'};

% xs and xw only depend on molality (not P or T); fnGval returns them as
% [1 × nm] in grid mode rather than broadcasting to [nP × nT × nm], so
% they cannot be assembled via assemble_grid.  Copy from either spline.
P_INDEP_FIELDS = {'xs', 'xw'};

Results = struct();

for k = 1:numel(all_fields)
    fn = all_fields{k};

    if any(strcmp(fn, COORD_FIELDS))
        continue;   % handled after the loop
    end

    if any(strcmp(fn, P_INDEP_FIELDS))
        if is_grid
            % Grid: xs/xw are [1 × nm] (P-independent) — copy from either spline
            if has_LP,  Results.(fn) = res_LP.(fn);
            else,       Results.(fn) = res_HP.(fn);   end
        else
            % Scatter: xs/xw are (N,) — must assemble across all N input points
            a_lp = [];  a_hp = [];
            if has_LP, a_lp = res_LP.(fn)(:); end
            if has_HP, a_hp = res_HP.(fn)(:); end
            Results.(fn) = assemble_scatter(a_lp, a_hp, nP, ...
                lp_only, lp_only_in, hp_only, hp_only_in, ...
                blend_global, blend_in_lp, blend_in_hp, w_blend);
        end
        continue;
    end

    if is_grid
        A_lp = [];  A_hp = [];
        if has_LP, A_lp = res_LP.(fn); end
        if has_HP, A_hp = res_HP.(fn); end
        Results.(fn) = assemble_grid(A_lp, A_hp, nP, ...
            lp_only, lp_only_in, hp_only, hp_only_in, ...
            blend_global, blend_in_lp, blend_in_hp, w_blend);
    else
        a_lp = [];  a_hp = [];
        if has_LP, a_lp = res_LP.(fn)(:); end
        if has_HP, a_hp = res_HP.(fn)(:); end
        Results.(fn) = assemble_scatter(a_lp, a_hp, nP, ...
            lp_only, lp_only_in, hp_only, hp_only_in, ...
            blend_global, blend_in_lp, blend_in_hp, w_blend);
    end
end

% ---- Echo coordinate fields directly from input -------------------------
if any(strcmp(all_fields, 'P'))
    if is_grid,  Results.P = P_vec(:);
    else,        Results.P = PT(:, 1);   end
end
if any(strcmp(all_fields, 'T'))
    if is_grid,  Results.T = PT{2}(:);
    else,        Results.T = PT(:, 2);   end
end
if any(strcmp(all_fields, 'm'))
    if is_grid,  Results.m = PT{3}(:);
    else,        Results.m = PT(:, 3);   end
end

end  % main function


% =========================================================================
%  Local helpers
% =========================================================================

function out = assemble_grid(A_lp, A_hp, nP, ...
                             lp_only, lp_only_in, hp_only, hp_only_in, ...
                             blend_idx, blend_in_lp, blend_in_hp, w_blend)
% Assemble a [nP × ...] array from LP and HP sub-arrays, blending where
% both splines overlap.  Trailing dimensions (T, m) are preserved.

% Use whichever result is available to infer trailing shape
if ~isempty(A_lp)
    sz = size(A_lp);
else
    sz = size(A_hp);
end
if numel(sz) < 2, sz = [sz, 1]; end
tail  = sz(2:end);
ntail = prod(tail);

out2 = zeros(nP, ntail);

if ~isempty(lp_only) && ~isempty(A_lp)
    lp2           = reshape(A_lp, size(A_lp, 1), ntail);
    out2(lp_only, :) = lp2(lp_only_in, :);
end
if ~isempty(hp_only) && ~isempty(A_hp)
    hp2           = reshape(A_hp, size(A_hp, 1), ntail);
    out2(hp_only, :) = hp2(hp_only_in, :);
end
if ~isempty(blend_idx) && ~isempty(A_lp) && ~isempty(A_hp)
    lp2 = reshape(A_lp, size(A_lp, 1), ntail);
    hp2 = reshape(A_hp, size(A_hp, 1), ntail);
    out2(blend_idx, :) = bsxfun(@times,   w_blend, lp2(blend_in_lp, :)) + ...
                         bsxfun(@times, 1-w_blend, hp2(blend_in_hp, :));
end

out = reshape(out2, [nP, tail]);
end


function out = assemble_scatter(a_lp, a_hp, nP, ...
                                lp_only, lp_only_in, hp_only, hp_only_in, ...
                                blend_idx, blend_in_lp, blend_in_hp, w_blend)
% Assemble a [nP × 1] vector from LP and HP sub-vectors.

out = zeros(nP, 1);

if ~isempty(lp_only) && ~isempty(a_lp)
    out(lp_only) = a_lp(lp_only_in);
end
if ~isempty(hp_only) && ~isempty(a_hp)
    out(hp_only) = a_hp(hp_only_in);
end
if ~isempty(blend_idx) && ~isempty(a_lp) && ~isempty(a_hp)
    out(blend_idx) = w_blend .* a_lp(blend_in_lp) + ...
                     (1 - w_blend) .* a_hp(blend_in_hp);
end
end
