function [colloc, b, ids] = collocate(knots, k, x, nd)
% COLLOCATE  Evaluate B-spline basis functions and their derivatives at data sites.
%
%   [colloc, b, ids] = collocate(knots, k, x, nd)
%
%   Inputs
%     knots  - knot sequence (length n+k) with multiplicity k at each end
%     k      - spline order (degree + 1)
%     x      - data sites (must be non-decreasing and within knot range)
%     nd     - number of derivative levels requested (0 = function only,
%              1 = function + 1st derivative, etc.)
%
%   Outputs
%     colloc - array of size (numDataPoints, numControls, nd+1)
%                colloc(:,:,1)  = basis functions (B-spline values)
%                colloc(:,:,i)  = (i-1)-th derivative,  i = 2 .. nd+1
%     b      - raw de Boor recurrence block, size (numDerivatives*numDataPoints, k)
%                Row layout: [pt1_deriv1; pt1_deriv2; ...; pt2_deriv1; ...]
%     ids    - (numDataPoints x k) array of 1-based column indices into colloc;
%                ids(j,:) gives the k active control-point columns for data point j
%
%   Algorithm: de Boor recurrence, vectorised simultaneously over all data points.
%   Reference: C. de Boor, "A Practical Guide to Splines", Algorithm 5.1 (Chapter X).
%
%   JMB 2019-2026.  No toolbox dependencies.
 
    numKnots    = length(knots);
    numControls = numKnots - k;       % = number of B-spline coefficients
    pts         = x(:);               % ensure column vector
    numDataPoints = length(pts);
    knots       = knots(:);
    numDerivatives = nd + 1;          % total levels: 1 = function only
    km1         = k - 1;
 
    % ------------------------------------------------------------------
    % Input validation
    % ------------------------------------------------------------------
    dk = diff(knots);
    idx_pos = find(dk > 0);
 
    if min(dk) < 0
        error('collocate:knotsNotNondecreasing', ...
              'Knot sequence must be non-decreasing.');
    end
    if numDerivatives > k
        error('collocate:derivativeTooHigh', ...
              'Requested derivative order (%d) exceeds spline order (%d).', nd, k);
    end
    if idx_pos(1) ~= k
        error('collocate:wrongStartMultiplicity', ...
              'Knot multiplicity at start must equal spline order k=%d.', k);
    end
    if idx_pos(end) ~= numControls
        error('collocate:wrongEndMultiplicity', ...
              'Knot multiplicity at end must equal spline order k=%d.', k);
    end
    if min(diff(x)) < 0
        error('collocate:xNotNondecreasing', ...
              'Data sites x must be non-decreasing.');
    end
    if x(1) < knots(1)
        error('collocate:xBelowKnotRange', ...
              'First data site lies below the first knot.');
    end
    if x(end) > knots(end)
        error('collocate:xAboveKnotRange', ...
              'Last data site lies above the last knot.');
    end
 
    % ------------------------------------------------------------------
    % Find the active knot span for each data point (savl).
    %
    % savl(j) = index i such that knots(i) <= x(j) < knots(i+1),
    % clamped to [k, numControls].  Uses a merge-sort trick: interleave
    % the mesh sites (interior knots) with pts, sort, then count how many
    % mesh sites precede each pt.
    % ------------------------------------------------------------------
    meshsites = knots(1:numControls);
    [~, sort_idx] = sort([meshsites(:).', pts(:).']);
    savl = max(find(sort_idx > numControls) - (1:numDataPoints), k);
 
    % ids(j, 1:k): 1-based column indices of the k active basis functions
    % for data point j.  Active functions span columns savl(j)-km1 .. savl(j).
    ids = savl(:) - km1 + (0:km1);   % (numDataPoints x k)
 
    % ------------------------------------------------------------------
    % de Boor recurrence to fill b
    %
    % b is a (numDerivatives*numDataPoints) x k block matrix.
    % Rows are interleaved by point then derivative:
    %   row  i + numDerivatives*(j-1)  corresponds to data point j,
    %   derivative level i  (1-based).
    % ------------------------------------------------------------------
    b = zeros(numDataPoints, k);
 
    if numDerivatives == 1
        % ---- No derivatives needed: standard B-spline recurrence ----
        % Initialise with the degree-0 basis (piecewise constant = 1 in
        % the active span).
        b(:, 1) = ones(numDataPoints, 1);
 
        for j = 1:km1
            saved = zeros(numDataPoints, 1);
            for r = 1:j
                tr   = knots(savl + r) - pts;       % right knot distance
                tl   = pts - knots(savl + r - j);   % left  knot distance
                term = b(:, r) ./ (tr + tl);
                b(:, r) = saved + tr .* term;
                saved   = tl .* term;
            end
            b(:, j + 1) = saved;
        end
 
    else
        % ---- Derivatives required: extended recurrence ----
        % bb interleaves numDerivatives rows per data point.
        % lptss tracks which rows of bb correspond to the "current" level.
        bb    = repmat([1, zeros(1, km1)], numDerivatives * numDataPoints, 1);
        lptss = numDerivatives * (1:numDataPoints);   % row indices of last deriv level
 
        % Phase 1: recurrence up to degree (k - numDerivatives), building
        % all derivative levels simultaneously.
        for j = 1:k - numDerivatives
            saved = zeros(numDataPoints, 1);
            for r = 1:j
                tr   = knots(savl + r) - pts;
                tl   = pts - knots(savl + r - j);
                term = bb(lptss, r) ./ (tr + tl);
                bb(lptss, r) = saved + tr .* term;
                saved        = tl .* term;
            end
            bb(lptss, j + 1) = saved;
        end
 
        % Phase 2: save derivative sub-blocks level by level.
        for jj = 1:numDerivatives - 1
            j     = k - numDerivatives + jj;
            saved = zeros(numDataPoints, 1);
            lptsn = lptss - 1;                        % row index for next level up
            for r = 1:j
                tr   = knots(savl + r) - pts;
                tl   = pts - knots(savl + r - j);
                term = bb(lptss, r) ./ (tr + tl);
                bb(lptsn, r) = saved + tr .* term;
                saved        = tl .* term;
            end
            bb(lptsn, j + 1) = saved;
            lptss = lptsn;
        end
 
        % Phase 3: convert saved B-spline blocks to derivative values by
        % differencing coefficients (de Boor, Chapter X, eq. 10.12).
        for jj = numDerivatives - 1:-1:1
            j       = k - jj;
            tempIdx = repmat((jj:numDerivatives-1)', 1, numDataPoints) ...
                    + repmat(lptsn, numDerivatives - jj, 1);
            lptss   = tempIdx(:);
            for r = j:-1:1
                tempDiff = repmat((knots(savl + r) - knots(savl + r - j)).' / j, ...
                                  numDerivatives - jj, 1);
                bb(lptss, r)   = -bb(lptss, r) ./ tempDiff(:);
                bb(lptss, r+1) =  bb(lptss, r+1) - bb(lptss, r);
            end
        end
        b = bb;
    end
 
    % ------------------------------------------------------------------
    % Scatter b into the sparse 3-D colloc array  (vectorised, no loops)
    %
    % b row layout:  row = deriv + numDerivatives*(pt-1)   (1-based)
    %
    % For each (deriv, pt) pair and each of the k basis columns c:
    %   colloc(pt, ids(pt,c), deriv) = b(deriv + numDerivatives*(pt-1), c)
    %
    % Build flat index arrays, then assign in one vectorised statement.
    % ------------------------------------------------------------------
    colloc = zeros(numDataPoints, numControls, numDerivatives);
 
    % Enumerate all (deriv, pt) combinations -- MATLAB column-major order
    [deriv_grid, pt_grid] = ndgrid(1:numDerivatives, 1:numDataPoints);
    b_rows = deriv_grid(:) + numDerivatives * (pt_grid(:) - 1); % row in b
 
    % Expand each (deriv, pt) pair over the k basis columns
    b_rows_exp  = repmat(b_rows,       1, k);   % (nD*nPts) x k
    deriv_exp   = repmat(deriv_grid(:),1, k);   % (nD*nPts) x k
    pt_exp      = repmat(pt_grid(:),   1, k);   % (nD*nPts) x k
    col_sel     = repmat(1:k, numel(b_rows), 1);% (nD*nPts) x k : basis-column index
 
    % ids(pt, col) gives the target control-point column in colloc
    ctrl_exp = ids(pt_exp + numDataPoints * (col_sel - 1));  % linear index into ids
 
    % Gather values from b using (row, col) linear index
    b_lin = b_rows_exp + size(b, 1) * (col_sel - 1);
    b_vals = b(b_lin);
 
    % Build linear index into colloc and scatter
    lin_colloc = pt_exp + numDataPoints * (ctrl_exp - 1) ...
               + numDataPoints * numControls * (deriv_exp - 1);
    colloc(lin_colloc(:)) = b_vals(:);
 
end
 