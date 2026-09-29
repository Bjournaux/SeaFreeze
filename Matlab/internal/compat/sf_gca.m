function ax = sf_gca(fig)
%SF_GCA  Current axes of a figure, creating one if the figure has none.
%
%   Portable stand-in for MATLAB's gca(fig); Octave's gca takes no argument.
%
%   See also: sf_ishandle_type

ax = get(fig, 'currentaxes');
if isempty(ax)
    ax = axes('Parent', fig);
end
end
