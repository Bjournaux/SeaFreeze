function tf = sf_ishandle_type(x, type)
%SF_ISHANDLE_TYPE  True when x is a live graphics handle of the given type.
%
%   Portable stand-in for isgraphics(x, type), whose availability varies
%   across Octave versions.  Like isgraphics, a bare number that happens to
%   match an existing handle counts as one.
%
%   See also: sf_gca

tf = false;
if isempty(x)
    return
end
try
    if ~all(ishghandle(x(:)))
        return
    end
    t = get(x, 'type');
    if ischar(t)
        tf = strcmp(t, type);
    else
        tf = all(strcmp(t, type));
    end
catch
    tf = false;
end
end
