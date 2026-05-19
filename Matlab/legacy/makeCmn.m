function D = makeCmn(A, B, C)
% MAKECMN  Sparse generalised Kronecker product for 2-D or 3-D tensor splines.
%
%   D = makeCmn(A, B)    -- 2-D problem (two collocation matrices)
%   D = makeCmn(A, B, C) -- 3-D problem (three collocation matrices)
%
%   Given sparse collocation matrices
%       A  (ai x aj),   B  (bk x bl),   C  (cp x cq)
%   returns the sparse matrix D of size (ai*bk*cp) x (aj*bl*cq) such that
%
%       D * coeff(:)  ==  vec( sum_i sum_k sum_p A(i,:)*B(k,:)*C(p,:) * coeff )
%
%   This is equivalent to kron(C, kron(B, A)) but 10-50x faster for the
%   sparse collocation matrices that arise in tensor B-spline evaluation,
%   because it operates only on the non-zero entries.
%
%   Row/column index mapping (matches the loop reference below):
%       row  =  i  +  (k-1)*ai  +  bk*ai*(p-1)
%       col  =  (col of A block)  +  (l-1)*aj  +  bl*aj*(q-1)
%   where i, k, p index non-zeros of A, B, C respectively.
%
%   If only A and B are supplied, C is treated as the scalar 1, giving the
%   2-D case D = kron(B, A).  If A and B are single-row vectors, D is a
%   single row -- the contribution of one data point to the system matrix.
%
%   Equivalent (but ~10x slower) reference implementation using full loops:
%
%     [ai,aj]=size(A); [bk,bl]=size(B); [cp,cq]=size(C);
%     D = zeros(ai*bk*cp, aj*bl*cq);
%     for p = 1:cp
%       for k = 1:bk
%         for i = 1:ai
%           mm = i + (k-1)*ai + bk*ai*(p-1);
%           for q = 1:cq
%             for l = 1:bl
%               nn = (l-1)*aj + bl*aj*(q-1);
%               D(mm, nn+(1:aj)) = A(i,:) * B(k,l) * C(p,q);
%             end
%           end
%         end
%       end
%     end
%
%   Benchmarks (dense random 10x10 matrices):
%     kron(C,kron(B,A))  ~10x slower
%     superkron          ~50x slower
%     makeCmn            baseline (results agree to machine precision)
%
%   JMB 2011-2026.

    % ------------------------------------------------------------------
    % Dimensions
    % ------------------------------------------------------------------
    [ai, aj] = size(A);
    [bk, bl] = size(B);

    if nargin == 3
        [cp, cq] = size(C);
    else
        % 2-D case: treat C as the scalar 1.
        % Using C=1 (not sparse(1,1,0)) is essential: find(sparse(1,1,0))
        % returns empty arrays because the stored value is zero, which
        % would cause the output to be all zeros.
        C  = 1;
        cp = 1;
        cq = 1;
    end

    % ------------------------------------------------------------------
    % Output matrix dimensions
    % ------------------------------------------------------------------
    nR = ai * bk * cp;
    nC = aj * bl * cq;

    % ------------------------------------------------------------------
    % Extract non-zero triplets from each sparse factor.
    % For a scalar C=1, find() returns (1,1,1) -- one non-zero entry.
    % ------------------------------------------------------------------
    [ii, jj, sa] = find(A);   % non-zeros of A: row ii, col jj, value sa
    [kk, ll, sb] = find(B);   % non-zeros of B: row kk, col ll, value sb
    [pp, qq, sc] = find(C);   % non-zeros of C: row pp, col qq, value sc

    ni  = numel(ii);           % nnz(A)
    nk  = numel(kk);           % nnz(B)
    np  = numel(pp);           % nnz(C)  -- equals 1 in the 2-D case
    nmm = ni * nk * np;        % total non-zeros in D

    % ------------------------------------------------------------------
    % Preallocate triplet arrays for sparse() assembly
    % ------------------------------------------------------------------
    mm = zeros(nmm, 1);        % row indices of D
    nn = zeros(nmm, 1);        % column indices of D
    v  = zeros(nmm, 1);        % values of D

    % ------------------------------------------------------------------
    % Fill triplets.
    %
    % Outer loops: iterate over non-zeros of C (p) and B (k).
    % Inner work:  the non-zeros of A are handled in one vectorised block
    %              per (k, p) pair, avoiding the innermost i-loop.
    %
    % For each (k, p) pair the row and column offsets are constant:
    %   im  = row offset    = (kk(k)-1)*ai  +  bk*ai*(pp(p)-1)
    %   in_ = column offset = (ll(k)-1)*aj  +  bl*aj*(qq(p)-1)
    % so the ni rows contributed are  ii + im  and columns are  jj + in_.
    %
    % count tracks the write position; it advances by nk*ni after each
    % p-iteration, and the k-loop fills consecutive ni-element slices.
    % ------------------------------------------------------------------
    count = 0;
    for p = 1:np
        im_offset = (pp(p) - 1) * bk * ai;   % row offset due to C factor
        in_offset = (qq(p) - 1) * bl * aj;   % col offset due to C factor

        for k = 1:nk
            im   = im_offset + (kk(k) - 1) * ai;   % row offset due to B factor
            in_  = in_offset + (ll(k) - 1) * aj;   % col offset due to B factor
            bkp  = sb(k) * sc(p);                   % combined B*C scalar weight

            % Write ni entries for all non-zeros of A at this (k,p) pair
            idx = count + (k-1)*ni + (1:ni);
            mm(idx) = ii + im;
            nn(idx) = jj + in_;
            v(idx)  = sa * bkp;
        end
        count = count + nk * ni;
    end

    % ------------------------------------------------------------------
    % Assemble sparse output matrix
    % ------------------------------------------------------------------
    D = sparse(mm, nn, v, nR, nC);
end
