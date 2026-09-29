function bin = sf_discretize(x, edges)
%SF_DISCRETIZE  Bin index of each element of x, as MATLAB's discretize.
%
%   bin = sf_discretize(x, edges) returns, for every element of x, the index i
%   with edges(i) <= x < edges(i+1).  The final bin is closed on the right, so
%   x == edges(end) maps to numel(edges)-1.  Values outside [edges(1),
%   edges(end)] and NaNs return NaN.
%
%   discretize is MATLAB-only.  histc has the binning behaviour we need and
%   exists in both interpreters (sp_val already relies on it for the gridded
%   evaluation path), so it is used here rather than a hand-rolled search.
%
%   edges must be non-decreasing.
%
%   See also: sf_is_octave

edges = edges(:).';
nbin  = numel(edges) - 1;
bin   = nan(size(x));
if nbin < 1
    return
end

% histc: edges(k) <= x < edges(k+1) -> k;  x == edges(end) -> numel(edges);
%        x out of range or NaN      -> 0.
[~, idx] = histc(x(:), edges);   %#ok<HISTC>
idx = reshape(idx, size(x));

inside      = (idx >= 1) & (idx <= nbin);
bin(inside) = idx(inside);
bin(idx == nbin + 1) = nbin;     % right edge closes the last bin
end
