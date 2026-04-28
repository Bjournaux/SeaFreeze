function out = SeaFreeze(varargin)
% SeaFreeze  (deprecated alias of SF_getprop, kept for backward compatibility)
% Baptiste Journaux - 2026
%
% As of v1.1.0 the primary entry point is SF_getprop. This wrapper forwards
% all arguments to SF_getprop and emits a one-time deprecation warning per
% MATLAB session. The warning can be suppressed with:
%
%     warning('off', 'SeaFreeze:deprecated')
%
% Plan to migrate calls of the form
%     out = SeaFreeze(PT, material, props)
% to
%     out = SF_getprop(PT, material, props)
%
% See `help SF_getprop` for full documentation.

persistent warned
if isempty(warned)
    warning('SeaFreeze:deprecated', ...
            ['SeaFreeze() is deprecated and will be removed in a future ' ...
             'release. Use SF_getprop() instead. ' ...
             '(Suppress with: warning(''off'',''SeaFreeze:deprecated''))']);
    warned = true;
end

out = SF_getprop(varargin{:});
end
