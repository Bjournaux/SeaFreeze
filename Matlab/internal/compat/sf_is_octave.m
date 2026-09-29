function tf = sf_is_octave()
%SF_IS_OCTAVE  True when running under GNU Octave rather than MATLAB.
%
%   Used by the handful of places in SeaFreeze where the two interpreters
%   genuinely differ (spline file format, missing builtins).  The answer
%   cannot change within a session, so it is cached.
%
%   See also: sf_ischarlike, sf_discretize, sf_load_spline

persistent cached
if isempty(cached)
    cached = (exist('OCTAVE_VERSION', 'builtin') == 5);
end
tf = cached;
end
