function [knots,k] = setknt(tau,k)
% SETKNT  Acceptable knot sequence satisfying Schoenberg-Whitney conditions.
%
%   SETKNT(TAU,K)  returns, for a given nondecreasing sequence TAU with
%   TAU(i) < TAU(i+K-1), all i, a knot sequence KNOTS for which the
%   Schoenberg-Whitney conditions
%
%            KNOTS(i) <  TAU(i)  <  KNOTS(i+K),   i = 1:length(TAU)
%
%   hold (with equality only at the first or last knot), ensuring that
%   the space of splines of order
%              K  :=  min(K, length(TAU))
%   with knot sequence KNOTS has a unique interpolant to arbitrary data
%   at the data sites TAU.  The K actually used is optionally returned.
%
%   This is the local replacement for the Spline Toolbox's APTKNT.
%   The knot sequence is also the initial guess used by OPTKNT.
%
%   See also AUGKNT, AVEKNT, NEWKNT, OPTKNT.
%
%   Algorithm: de Boor, "A Practical Guide to Splines" (2001), Ch. XII.
%   JMB 2019-2026  (based on MathWorks aptknt, Copyright 1987-2008).

n = length(tau);
if n < 2
    error('setknt:tooFewPoints', ...
          'TAU must have at least 2 points; got %d.', n);
end

k = max(1, min(k, n));
dtau = diff(tau);

if any(dtau < 0)
    error('setknt:tauDecreasing', ...
          'TAU must be non-decreasing.');
end

if k == 1
    % Order 1: use midpoints between data sites as knots.
    if ~all(dtau)
        error('setknt:tauNotStrictlyIncreasing', ...
              'TAU must be strictly increasing when k=1.');
    end
    knots = [tau(1), tau(1:n-1) + dtau/2, tau(n)];
else
    if any(tau(k:n) == tau(1:n-k+1))
        error('setknt:tauMultiplicityTooLarge', ...
              'TAU has too many repeated values for order k=%d (max %d consecutive equal values).', ...
              k, k-1);
    end
    knots = addknt([tau(1), avknt(tau,k), tau(end)], k);
end



function [augknot,addl] = addknt(knots,k,mults)
% ADDKNT  Augment a knot sequence to have multiplicity k at each end.
%
%   ADDKNT(KNOTS,K) returns a nondecreasing augmented knot sequence
%   with the first and last knot repeated exactly K times.
%
%   ADDKNT(KNOTS,K,MULTS) also repeats each interior knot MULTS times.
%   MULTS may be a scalar (uniform) or a vector (one entry per interior knot).
%
%   Private helper for setknt.  Equivalent to the Toolbox's AUGKNT.

if nargin < 3
    if (length(k) > 1 || k < 1)   % || not |: short-circuit for scalar guard
        error('setknt:addknt:badK', 'k must be a positive scalar integer.');
    end
    mults = 1;
end

dk = diff(knots);
if ~isempty(find(dk < 0, 1))
    knots = sort(knots);
    dk    = diff(knots);
end

j = find(dk > 0);
if isempty(j)
    error('setknt:addknt:tooFewKnots', ...
          'Knot sequence has no positive-length intervals.');
end
addl = k - j(1);

interior = (j(1)+1) : j(end);
if length(mults) ~= length(interior)
    mults = repmat(mults(1), size(interior));
end

augknot = bk2kt(knots([1, interior, end]), [k, mults, k]);

function tstar = avknt(t,k)
% AVKNT  Knot averages: tstar(i) = mean(t(i+1:i+k-1)).
%
%   These are the recommended interpolation sites when interpolating
%   from the spline space S_{k,t}.  Private helper for setknt.

t = t(:);
n = length(t) - k;
if k < 2
    error('setknt:avknt:badK', 'k must be at least 2.');
elseif n < 0
    error('setknt:avknt:tooFewKnots', 'Too few knots for the requested order.');
elseif k == 2
    tstar = reshape(t(1+(1:n)), 1, n);
else
    temp  = repmat(t, 1, k-1);
    temp  = sum(reshape([temp(:); zeros(k-1,1)], n+k+1, k-1).')/(k-1);
    tstar = temp(1+(1:n));
end

function t = bk2kt(breaks,mults)
% BK2KT  Expand breaks with multiplicities into a full knot sequence.
%
%   T = BK2KT(BREAKS, MULTS) returns the sequence T in which BREAKS(i)
%   is repeated MULTS(i) times.  MULTS may be a scalar (uniform repeat
%   count) or a vector with one entry per break.
%
%   Example:
%      t = bk2kt([1 2], 3)   % gives [1 1 1 2 2 2]
%
%   Private helper for setknt/addknt.  Equivalent to the Toolbox's BRK2KNT.

s = sum(mults);
if s==0
    t = [];
else
    li = length(breaks);
    % make sure there is a multiplicity assigned to each break,
    % and drop any break whose assigned multiplicity is not positive.
    if length(mults)~=li 
        mults = repmat(mults(1),1,li);
        s = mults(1)*li;
    else
        fm = find(mults<=0);
        if ~isempty(fm)
            breaks(fm)=[];
            mults(fm)=[];
            li = length(breaks);
        end
    end
    mm = zeros(1,s);
    mm(cumsum([1 reshape(mults(1:li-1),1,li-1)])) = ones(1,li);
    t = breaks(cumsum(mm));
end

