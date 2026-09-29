function tf = sf_istextlike(x)
%SF_ISTEXTLIKE  True for a char array or a string array of any size.
%
%   Portable stand-in for the MATLAB idiom
%       ischar(x) || isstring(x)
%   used where a whole string array is acceptable (e.g. a list of property
%   names that is about to be passed through cellstr).  Octave has no string
%   class, so only the char test applies there.
%
%   See also: sf_ischarlike, sf_is_octave

tf = ischar(x);
if ~tf && ~sf_is_octave()
    tf = isstring(x);
end
end
