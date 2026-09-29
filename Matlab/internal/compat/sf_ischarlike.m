function tf = sf_ischarlike(x)
%SF_ISCHARLIKE  True for a character vector or a scalar string.
%
%   Portable stand-in for the MATLAB idiom
%       ischar(x) || (isstring(x) && isscalar(x))
%   Octave has no string class and no isstring, so the second test is only
%   attempted under MATLAB.
%
%   See also: sf_istextlike, sf_is_octave

tf = ischar(x);
if ~tf && ~sf_is_octave()
    tf = isstring(x) && isscalar(x);
end
end
