function cmap = sf_parula(n)
%SF_PARULA  MATLAB's parula colormap, or the closest equivalent on Octave.
%
%   Octave does not implement parula. It does ship viridis, which is the same
%   kind of perceptually uniform dark-blue-to-yellow ramp, so that is used
%   instead: Octave figures will not be pixel-identical to MATLAB's, but the
%   curves stay distinguishable and ordered the same way. jet is a last resort.
%
%   See also: sf_is_octave

if nargin < 1
    n = 64;
end

if exist('parula', 'file') || exist('parula', 'builtin')
    cmap = parula(n);
elseif exist('viridis', 'file') || exist('viridis', 'builtin')
    cmap = viridis(n);
else
    cmap = jet(n);
end
end
