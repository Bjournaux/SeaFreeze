function out = sp_val(sp, derv, x)
% SP_VAL  Evaluate a B-form spline (tensor or univariate) without the Toolbox.
%
%   out = sp_val(sp, derv, x)
%   out = sp_val(sp, x)          % shorthand when no derivatives needed
%
%   Substitutes for fnval + fnder from the MATLAB Spline Toolbox.
%   Self-contained: calls no external or toolbox functions.
%
%   Inputs
%     sp    B-form spline struct (fields: knots, coefs, number, order, dim)
%     derv  Derivative order vector, one element per independent variable.
%           sp_val(sp, [0 1 0], x) gives d/dT for a (P,T,m) spline.
%           Omit or pass zeros for function values only.
%     x     Evaluation points:
%             Gridded:   cell array {x1, x2, ...} -- each xi is a vector.
%             Scattered: matrix of size (nPoints x nDims).
%
%   Output
%     out   Evaluated values.
%             Gridded:   array of size (n1 x n2 x ...) or (d x n1 x n2 x ...)
%             Scattered: column vector of length nPoints (or d x nPoints if d>1)
%
%   Algorithm
%     Gridded data: iterative B-form evaluation along each dimension.
%     Scattered data: B-form converted to pp-form via sp2pp, then evaluated
%     at all points simultaneously via fully-vectorised ppvl.  The pp struct
%     is cached (persistent) so repeated calls with the same differentiated
%     spline (e.g. multiple derivatives from fnGval) skip the conversion.
%
%   References
%     C. de Boor, "A Practical Guide to Splines", Springer, 2001.
%
%   JMB 2019-2026.

if nargin == 2
    x    = derv;
    derv = zeros(1, size(x, 2));    % default: no derivatives
end

sp  = fndr(sp, derv);    % differentiate by coefficient differencing
out = spvl(sp, x);
end


% ======================================================================
%  EVALUATION DISPATCHER
% ======================================================================

function out = spvl(sp, x)
% SPVL  Dispatch B-form evaluation to gridded or scattered path.

% Persistent pp cache for scattered evaluation.
% Key = hash of sp.coefs (using a fast checksum).
% This avoids re-running sp2pp for every derivative call on the same spline.
% Cache is cleared automatically when sp_val.m is edited or 'clear functions'.
persistent pp_store

if isempty(pp_store)
    pp_store = struct('key', {}, 'pp', {});
end

if iscell(sp.knots)
    % ---- Tensor-product spline ----
    [t, a, n, ~, d] = spbk(sp);
    m = length(t);

    if iscell(x)
        % ---- Gridded evaluation: iterative B-form along each dimension ----
        % Efficient for large regular grids; no conversion needed.
        v = a;
        sizev    = [d, n];
        nsizev   = zeros(1, m);
        for i = m:-1:1
            nsizev(i) = length(x{i}(:));
            v = reshape(spv(spmk(t{i}, reshape(v, prod(sizev(1:m)), sizev(m+1))), ...
                            x{i}), [sizev(1:m), nsizev(i)]);
            sizev(m+1) = nsizev(i);
            if m > 1
                v = permute(v, [1, m+1, 2:m]);
                sizev(2:m+1) = sizev([m+1, 2:m]);
            end
        end
        if d > 1
            out = reshape(v, [d, nsizev]);
        else
            out = reshape(v, nsizev);
        end

    else
        % ---- Scattered evaluation: convert to pp-form then evaluate ----
        % Look up or compute the pp struct for this (differentiated) spline.
        cache_key = sp_cache_key(sp);
        pp = [];
        for ci = 1:numel(pp_store)
            if strcmp(pp_store(ci).key, cache_key)
                pp = pp_store(ci).pp;
                break;
            end
        end
        if isempty(pp)
            pp = sp2pp(sp);
            pp_store(end+1) = struct('key', cache_key, 'pp', pp);
            % Keep cache bounded: evict oldest entry beyond 32 entries
            if numel(pp_store) > 32
                pp_store = pp_store(2:end);
            end
        end
        out = ppvl(pp, x);
    end

else
    % ---- Univariate spline ----
    if iscell(x)
        out = spv(sp, x{1});
    else
        out = spv(sp, x);
    end
end
end


function key = sp_cache_key(sp)
% Build a compact string key from the spline coef content.
% Uses a fast sum+xor checksum -- not cryptographic but sufficient
% for distinguishing different differentiated versions of a spline.
c = sp.coefs(:);
n = numel(c);
% Combine: sum, max, and a strided sample for speed on large arrays
stride = max(1, floor(n/200));
sample = c(1:stride:end);
key = sprintf('%.6g_%.6g_%.6g_%d', sum(sample), max(abs(sample)), ...
              sample(ceil(end/2)), n);
end


% ======================================================================
%  B-FORM -> PP-FORM CONVERSION  (tensor and univariate)
% ======================================================================

function pp = sp2pp(sp)
% SP2PP  Convert a B-form spline to piecewise polynomial (pp) form.
%
%   Works for both univariate and tensor-product splines.
%   For tensor splines the conversion is applied dimension by dimension,
%   mirroring the gridded evaluation loop in spvl.
%
%   The resulting pp struct is compatible with ppvl and ppmk.

if iscell(sp.knots)
    % ---- Tensor spline: convert one dimension at a time ----
    [t, coefs, n, ~, d] = spbk(sp);
    m      = length(t);
    sizec  = [d, n];
    breaks = cell(1, m);

    % sp.coefs may be stored as a flat 2D matrix (d*n1 x n2) or as an ND
    % array (d x n1 x n2 x ...).  Canonicalise to the ND form [d, n1, ..., nm]
    % so that all subsequent reshapes are consistent.
    coefs = reshape(coefs, sizec);

    for i = m:-1:1
        % Reshape coefficients so the current dimension is last.
        % Mirrors the MathWorks sp2pp exactly.
        d_eff = prod(sizec(1:m));
        uni   = spmk(t{i}, reshape(coefs, d_eff, sizec(m+1)));
        ppi   = s2p1(uni);
        breaks{i}  = ppi.breaks;
        sizec(m+1) = ppi.pieces * ppi.order;
        coefs      = reshape(ppi.coefs, sizec);   % no extra permute of ppi.coefs
        if m > 1
            coefs        = permute(coefs, [1, m+1, 2:m]);
            sizec(2:m+1) = sizec([m+1, 2:m]);
        end
    end

    pp = ppmk(breaks, coefs, sizec);

else
    % ---- Univariate spline ----
    pp = s2p1(sp);
end
end


% ======================================================================
%  PP-FORM EVALUATOR AT SCATTERED POINTS
% ======================================================================

function out = ppvl(pp, x)
% PPVL  Evaluate a pp-form spline at scattered points.
%
%   pp   pp struct from ppmk (univariate or tensor).
%   x    Scattered points:
%          univariate: column vector or row vector
%          tensor:     matrix (nPoints x nDims)
%
%   Evaluation uses Horner's method for each polynomial piece,
%   vectorised over all points simultaneously.
%
%   The pp coef convention (from ppmk / iUnivariateWithoutD) is:
%     coefs is (d*l) x k, where row d*(p-1)+comp holds the coefficients
%     for piece p and component comp, from highest to lowest power.
%   For tensor splines the coefs are stored as [d, l1*k1, l2*k2, ...].

if iscell(pp.breaks)
    out = ppvl_tensor(pp, x);
else
    out = ppvl_uni(pp, x(:));
end
end


function out = ppvl_uni(pp, x)
% PPVL_UNI  Evaluate a univariate pp spline at scattered x (column vector).
%
%   Finds the piece each point belongs to, then evaluates the local
%   polynomial via Horner's method, all vectorised over nPoints.

breaks = pp.breaks(:).';    % row vector, length l+1
k      = pp.order;
l      = pp.pieces;
d      = prod(pp.dim);
npts   = numel(x);

% Find piece index for each point: piece p means breaks(p) <= x < breaks(p+1)
% Clamp to [1, l] so boundary points land in the last piece.
p = min(discretize(x, breaks), l);

% Handle points at or beyond the right boundary
p(isnan(p) | p < 1) = 1;    % left of domain -> piece 1 (will zero out below)

% breaks is a row vector; indexing with column p(:) returns a row in MATLAB,
% so reshape to column explicitly.
dx = x(:) - reshape(breaks(p(:)), [], 1);   % (npts x 1)

% Extract coefficient rows for each point: row index = d*(p-1) + comp
% For scalar (d=1): row = p.  Generalise for d>1 below.
% coefs is (d*l) x k.
row_base = d * (p(:) - 1);    % (npts x 1), 0-based start for each piece

out = zeros(d, npts);
for comp = 1:d
    rows = row_base + comp;               % row indices into coefs (1-based)
    c    = pp.coefs(rows, :);             % (npts x k), highest power first
    % Horner evaluation: c(:,1)*dx^(k-1) + ... + c(:,k-1)*dx + c(:,k)
    val  = c(:, 1);
    for j = 2:k
        val = val .* dx + c(:, j);
    end
    out(comp, :) = val.';
end

% Zero out points outside the domain
outside = (x(:) < breaks(1)) | (x(:) > breaks(end));
out(:, outside) = 0;

if d == 1
    out = out(:);    % return column vector for scalar splines
end
end


function out = ppvl_tensor(pp, x)
% PPVL_TENSOR  Evaluate a tensor pp spline at scattered points.
%
%   x is (nPoints x m), one column per independent variable.
%
%   Coef layout (from MathWorks-compatible sp2pp): pp.coefs has shape
%   [d, l1*k1, l2*k2, ..., lm*km].  Within each li*ki block the ordering
%   is PIECE-MINOR (power-major):
%       column lk = piece + (power-1)*li   (1-based)
%
%   Algorithm: reduce coefs one dimension at a time (i = m down to 1).
%   C is maintained as (nc x npts_cur):
%     - First iteration: npts_cur=1, Cmat is 2-D (nc_next x li*ki).
%       Column gather is a simple matrix index: Cmat(:, cols(:,pw).')
%       giving (nc_next x npts) directly -- no loops.
%     - Subsequent iterations: npts_cur=npts, Cmat is 3-D.
%       Linear indexing gathers the correct element for each (row, pt)
%       pair without any explicit row or point loops.
%   Both steps use vectorised Horner evaluation.

breaks   = pp.breaks;
k        = pp.order;          % [k1, k2, ..., km]
l        = pp.pieces;         % [l1, l2, ..., lm]
m        = length(breaks);
d        = prod(pp.dim);
npts     = size(x, 1);

nc_total = d * prod(l .* k);
C        = reshape(pp.coefs, nc_total, 1);   % (nc_total x 1)
npts_cur = 1;

for i = m:-1:1
    li = l(i);  ki = k(i);
    br = breaks{i}(:).';
    xi = x(:, i);

    % Piece index (1-based, clamped) and local coordinate
    p  = min(discretize(xi, br), li);
    p(isnan(p) | p < 1) = 1;
    dx = xi(:) - reshape(br(p(:)), [], 1);   % (npts x 1)

    % Piece-minor column indices for each point and each power:
    %   cols(pt, pw) = p(pt) + (pw-1)*li
    cols    = p(:) + (0:ki-1)*li;             % (npts x ki)
    nc_next = round(numel(C) / (li * ki * npts_cur));
    Cmat    = reshape(C, nc_next, li*ki, npts_cur);

    if npts_cur == 1
        % ---- First iteration: Cmat is 2-D (nc_next x li*ki) ----
        % Cmat(:, v) where v is (npts x 1) returns (nc_next x npts).
        % Vectorised Horner: no row or point loops.
        result = Cmat(:, cols(:,1).');           % gather power 1
        for pw = 2:ki
            result = result .* dx(:).' + Cmat(:, cols(:,pw).');
        end

    else
        % ---- Subsequent iterations: Cmat is 3-D (nc_next x li*ki x npts) ----
        % Build linear index: lin(row,pt) = row + (col-1)*nc_next + (pt-1)*nc_next*li*ki
        % where col = cols(pt, pw).  Precompute the pt-dependent offset once.
        row_idx = repmat((1:nc_next).', 1, npts);          % (nc_next x npts)
        pt_off  = (0:npts-1) * (nc_next * li * ki);        % (1 x npts)
        pt_idx  = repmat(pt_off, nc_next, 1);              % (nc_next x npts)

        col_idx = repmat(cols(:,1).', nc_next, 1);         % (nc_next x npts)
        result  = Cmat(row_idx + (col_idx-1)*nc_next + pt_idx);
        for pw = 2:ki
            col_idx = repmat(cols(:,pw).', nc_next, 1);
            result  = result .* dx(:).' + ...
                      Cmat(row_idx + (col_idx-1)*nc_next + pt_idx);
        end
    end

    C        = result;       % (nc_next x npts)
    npts_cur = npts;
end

out = reshape(C, d, npts);

% Zero out points outside the valid domain in any dimension
outside = false(npts, 1);
for i = 1:m
    br      = breaks{i}(:).';
    outside = outside | (x(:,i) < br(1)) | (x(:,i) > br(end));
end
out(:, outside) = 0;
if d == 1, out = out(:); end
end


% ======================================================================
%  UNIVARIATE B-FORM EVALUATOR  (used for gridded path)
% ======================================================================

function v = spv(sp, x)
% SPV  Evaluate a univariate B-form spline at points x.
%   Based on de Boor Algorithm A (SPVAL1).
%   Used for gridded evaluation only; scattered data goes through ppvl.

[mx, nx] = size(x);
lx = mx * nx;
xs = reshape(x, 1, lx);

[t, a, n, k, d] = spbk(sp);
if lx == 0, v = zeros(d, 0); return, end

% Augment knot sequence so first and last knot have multiplicity >= k
index = find(diff(t) > 0);
addl  = k - index(1);
addr  = index(end) - n;
if addl > 0 || addr > 0
    npk = n + k;
    t   = t([ones(1,addl), 1:npk, npk*ones(1,addr)]);
    a   = [zeros(d, addl), a, zeros(d, addr)];
    n   = n + addl + addr;
end

% Find knot interval for each evaluation point using discretize
% (replaces deprecated histc)
xi_neg = -xs;
t_neg  = [-inf, -fliplr(t(k+1:n)), inf];
[~, index] = histc(xi_neg, t_neg);    % keep histc here: internal gridded path, not user-facing
NaNx  = find(index == 0);
index = max(n + 1 - index, k);
if ~isempty(NaNx), index(NaNx) = k; end

if k > 1
    % de Boor triangular algorithm, vectorised over all points
    dindex = reshape(repmat(index, d, 1), d*lx, 1);
    tx = reshape(t(repmat(2-k:k-1, d*lx, 1) + repmat(dindex, 1, 2*(k-1))), d*lx, 2*(k-1));
    tx = tx - repmat(reshape(repmat(xs, d, 1), d*lx, 1), 1, 2*(k-1));
    dindex = reshape(repmat(d*index, d, 1) + repmat((1-d:0).', 1, lx), d*lx, 1);
    b = repmat(d*(1-k):d:0, d*lx, 1) + repmat(dindex, 1, k);
    a = a(:); b(:) = a(b);
    for r = 1:k-1
        for i = 1:k-r
            b(:,i) = (tx(:,i+k-1) .* b(:,i) - tx(:,i+r-1) .* b(:,i+1)) ./ ...
                     (tx(:,i+k-1)            - tx(:,i+r-1));
        end
    end
    v = reshape(b(:,1), d, lx);
else
    v = a(:, index);
    if ~isempty(NaNx), v(:, NaNx) = NaN; end
end

% Zero out points outside the basic interval
outside = find(x < t(1) | x > t(n+k));
if ~isempty(outside)
    v(:, outside) = zeros(d, length(outside));
end
v = reshape(v, d*mx, nx);
end


% ======================================================================
%  SPLINE STRUCT HELPERS
% ======================================================================

function spline = spmk(knots, coefs, sizec)
% SPMK  Assemble a B-form spline struct.  Based on SPMAK.
if nargin < 3
    sizec = size(coefs);
end
m = 1;
if iscell(knots)
    m = length(knots);
end
if length(sizec) == m
    sizec = [1, sizec];
end
sizeval = sizec(1:end-m);
sizec   = [prod(sizeval), sizec(end-m+(1:m))];
coefs   = reshape(coefs, sizec);
if iscell(knots)
    [knots, coefs, k, sizec] = chkt(knots, coefs, sizec);
else
    [knots, coefs, k, sizec] = chkt({knots}, coefs, sizec);
    knots = knots{1};
end
spline.form   = 'B-';
spline.knots  = knots;
spline.coefs  = coefs;
spline.number = sizec(2:end);
spline.order  = k;
spline.dim    = sizeval;
end


function [knots, coefs, k, sizec] = chkt(knots, coefs, sizec)
% CHKT  Check knots and remove trivial B-splines.  Based on CHCKKNT.
for j = 1:length(sizec)-1
    n    = sizec(j+1);
    k(j) = length(knots{j}) - n;
    knots{j} = reshape(knots{j}, 1, n+k(j));
    index = find(knots{j}(k(j)+(1:n)) - knots{j}(1:n) > 0);
    if length(index) < n
        oldn = n; n = length(index);
        knots{j} = reshape(knots{j}([index, oldn+(1:k(j))]), 1, n+k(j));
        coefs = reshape(coefs, [prod(sizec(1:j)), sizec(j+1), prod(sizec(j+2:end))]);
        sizec(j+1) = n;
        coefs = reshape(coefs(:, index, :), sizec);
    end
end
end


function varargout = spbk(sp)
% SPBK  Break apart a spline struct into {t, a, n, k, d}.
varargout = {sp.knots, sp.coefs, sp.number, sp.order, sp.dim};
end


function fprime = fndr(f, dorder)
% FNDR  Differentiate a B-form spline by differencing coefficients.
%   Based on FNDER (de Boor).  Works for tensor and univariate splines.
sizeval = f.dim;
if length(sizeval) > 1, f.dim = prod(sizeval); end
if nargin < 2, dorder = 1; end

[knots, coefs, n, ~, d] = spbk(f);
if iscell(knots)
    m     = length(knots);
    sizec = [d, n];
    for i = m:-1:1
        dsp      = fndb(spmk(knots{i}, reshape(coefs, prod(sizec(1:m)), sizec(m+1))), dorder(i));
        knots{i} = dsp.knots;
        sizec(m+1) = dsp.number;
        coefs    = reshape(dsp.coefs, sizec);
        if m > 1
            coefs        = permute(coefs, [1, m+1, 2:m]);
            sizec(2:m+1) = sizec([m+1, 2:m]);
        end
    end
    fprime = spmk(knots, coefs, sizec);
else
    fprime = fndb(f, dorder);
end
if length(sizeval) > 1, fprime.dim = sizeval; end
end


function fprime = fndb(f, dorder)
% FNDB  Differentiate a univariate B-form spline by coefficient differencing.
%   Based on FNDERB (de Boor, Chapter X).
[t, a, n, k, d] = spbk(f);
if k <= dorder
    fprime = spmk(t, zeros(d, n));
elseif dorder < 0
    error('sp_val:fndb:integrationNotImplemented', ...
          'Integration (negative dorder) is not implemented.');
else
    knew = k - dorder;
    for j = k-1:-1:knew
        tt   = t(j+1+(0:n)) - t(1:n+1);
        z    = find(tt > 0);
        nn   = length(z);
        temp = diff([zeros(1,d); a.'; zeros(1,d)]).';
        a    = temp(:,z) ./ repmat(tt(z)/j, d, 1);
        t    = [t(z), t(n+2:n+j+1)];
        n    = nn;
    end
    fprime = spmk(t, a);
end
end


% ======================================================================
%  B-FORM -> PP-FORM (UNIVARIATE)
% ======================================================================

function pp = s2p1(spline)
% S2P1  Convert a univariate B-form spline to piecewise polynomial form.
%   Based on de Boor's FN2FM (B- to pp-form).
%   Uses knot insertion (srp) to obtain Taylor coefficients at each breakpoint.

[t, a, n, k, d] = spbk(spline);

% Augment knot sequence to multiplicity k at each end
index = find(diff(t) > 0);
addl  = k - index(1);
addr  = index(end) - n;
if addl > 0 || addr > 0
    t = [repmat(t(1), 1, addl), t(:).', repmat(t(n+k), 1, addr)];
    a = [zeros(d, addl), a, zeros(d, addr)];
end

% Identify interior breakpoints (knot intervals with positive length)
inter = find(diff(t) > 0);
l     = length(inter);

if k > 1
    % For each breakpoint, extract the local B-coefficients and convert
    % to Taylor coefficients via the knot-insertion algorithm srp.
    temp   = repmat(inter, d, 1);
    dinter = temp(:);
    tx     = repmat(2-k:k-1, d*l, 1) + repmat(dinter, 1, 2*(k-1));
    tx(:)  = t(tx);
    tx     = tx - repmat(t(dinter).', 1, 2*(k-1));
    a      = a(:);
    temp   = repmat(d*inter, d, 1) + repmat((1-d:0).', 1, l);
    dinter(:) = temp(:);
    b      = repmat(d*(1-k:0), d*l, 1) + repmat(dinter, 1, k);
    b(:)   = a(b);
    c      = srp(tx, b);
else
    c = a(:, inter);
    c = c(:);
end

pp = ppmk([t(inter), t(inter(end)+1)], c, d);
end


function [v, b] = srp(tx, a)
% SRP  Right Taylor coefficients from local B-coefficients via knot insertion.
%   Based on de Boor's SRP (25 Feb 1989).
%   See also: de Boor, "A Practical Guide to Splines", Algorithm A.
%
%   tx  (d*l x 2*(k-1)) knot differences relative to each breakpoint
%   a   (d*l x k)       local B-spline coefficients
%   v   (d*l x k)       Taylor coefficients, highest power last (reversed on output)
%   b   (d*l x k)       intermediate B-coefficients after knot insertion

k   = size(a, 2);
km1 = k - 1;
b   = a;

% Repeated knot insertion at 0 to get B-coefficients over [0, tx(:,k)]
for r = 1:km1
    for i = 1:k-r
        b(:,i) = (tx(:,i+km1) .* b(:,i) - tx(:,i+r-1) .* b(:,i+1)) ./ ...
                 (tx(:,i+km1)            - tx(:,i+r-1));
    end
end

% Differentiate at 0 to get Taylor (pp) coefficients
v = b;
for r = 2:k
    factor = (k - r + 1) / (r - 1);
    for i = k:-1:r
        v(:,i) = (v(:,i) - v(:,i-1)) * factor ./ tx(:,i+k-r);
    end
end

v = v(:, k:-1:1);    % reverse so highest power is first (ppmk convention)
end


% ======================================================================
%  PP-FORM STRUCT ASSEMBLER
% ======================================================================

function pp = ppmk(breaks, coefs, d)
% PPMK  Assemble a piecewise polynomial struct.
%   Based on MathWorks PPMAK.  Handles univariate and tensor-product splines.
%
%   For univariate:  coefs is (d x k*l); internally reshaped to (d*l x k).
%   For tensor:      coefs is [d, l1*k1, l2*k2, ...] (from sp2pp loop).

if nargin == 0
    breaks = input('Give the (l+1)-vector of breaks  > ');
    coefs  = input('Give the (d by (k*l)) matrix of coefficients  > ');
end

sizec = size(coefs);

if iscell(breaks)
    if nargin > 2
        if prod(sizec) ~= prod(d)
            error('sp_val:ppmk:coefsDontMatchSize', ...
                  'Size of coefs does not match the specified size d.');
        end
        sizec = d;
    end
    [breaks, coefs, sizeval, l, k] = iMultivariateSpline(breaks, coefs, sizec);
else
    if nargin < 3
        [coefs, sizeval, l, k] = iUnivariateWithoutD(breaks, coefs, sizec);
    else
        [coefs, sizeval, l, k] = iUnivariateWithD(breaks, coefs, sizec, d);
    end
    breaks = reshape(breaks, 1, l+1);
end

pp.form   = 'pp';
pp.breaks = breaks;
pp.coefs  = coefs;
pp.pieces = l;
pp.order  = k;
pp.dim    = sizeval;
end


function [breaks, coefs, sizeval, l, k] = iMultivariateSpline(breaks, coefs, sizec)
m = length(breaks);
if length(sizec) < m
    error('sp_val:ppmk:coefsLengthLessThanBreaks', ...
          'Length of sizec must be at least the number of break sequences.');
end
if length(sizec) == m
    sizec = [1, sizec];
end
sizeval = sizec(1:end-m);
sizec   = [prod(sizeval), sizec(end-m+(1:m))];
coefs   = reshape(coefs, sizec);
l = zeros(1, m); k = zeros(1, m);
for i = m:-1:1
    l(i) = length(breaks{i}) - 1;
    k(i) = fix(sizec(i+1) / l(i));
    if k(i) <= 0 || k(i)*l(i) ~= sizec(i+1)
        error('sp_val:ppmk:piecesDoNotMatchCoefs', ...
              'Dimension %d: %d pieces do not match coef size %d.', i, l(i), sizec(i+1));
    end
    breaks{i} = reshape(breaks{i}, 1, l(i)+1);
end
end


function [coefs, sizeval, l, k] = iUnivariateWithoutD(breaks, coefs, sizec)
if isempty(coefs)
    error('sp_val:ppmk:emptyCoefs', 'Coefficient array is empty.');
end
sizeval = sizec(1:end-1);
d       = prod(sizeval);
kl      = sizec(end);
l       = length(breaks) - 1;
k       = fix(kl / l);
if k <= 0 || k*l ~= kl
    error('sp_val:ppmk:piecesDoNotMatchCoefs', ...
          '%d pieces do not match coef length %d.', l, kl);
elseif any(diff(breaks) < 0)
    error('sp_val:ppmk:decreasingBreaks', 'Breaks must be non-decreasing.');
elseif breaks(1) == breaks(l+1)
    error('sp_val:ppmk:extremeBreaksSame', 'First and last break must differ.');
end
coefs = reshape(permute(reshape(coefs, [d, k, l]), [1, 3, 2]), d*l, k);
end


function [coefs, sizeval, l, k] = iUnivariateWithD(breaks, coefs, sizec, d)
if length(d) == 1
    k = sizec(end);
    l = prod(sizec(1:end-1)) / d;
else
    if prod(d) ~= prod(sizec)
        error('sp_val:ppmk:coefsSizeMismatch', ...
              'coefs size %s does not match specified size %s.', ...
              mat2str(sizec), mat2str(d));
    end
    k    = d(end);
    l    = d(end-1);
    d(end-1:end) = [];
    if isempty(d), d = 1; end
    coefs = reshape(coefs, prod(d)*l, k);
end
if l + 1 ~= length(breaks)
    error('sp_val:ppmk:coefsDontMatchBreaks', ...
          '%d pieces but %d break intervals.', l, length(breaks)-1);
end
sizeval = d;
end
