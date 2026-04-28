function sf_validate_PT(PT, material)
% Validate the PT argument shape and contents for SeaFreeze / SF_WhichPhase.
% Baptiste Journaux - 2026
% Throws 'SeaFreeze:badInput' on failure with a clear message.
%
%   PT       - cell {P,T} / {P,T,m} (gridded) or numeric N-by-2 / N-by-3 (scatter)
%   material - phase name; 'NaClaq' triggers the 3D (P,T,m) requirement,
%              everything else expects 2D (P,T).

is_solution = strcmp(material, 'NaClaq');
expected_dim = 2 + double(is_solution);
shape_str = ternary(is_solution, '{P,T,m} (cell) or N-by-3 (matrix)', ...
                                 '{P,T} (cell) or N-by-2 (matrix)');

if iscell(PT)
    if ~isvector(PT) || length(PT) ~= expected_dim
        error('SeaFreeze:badInput', ...
              ['Cell PT for material ''%s'' must have %d entries; got %d. ' ...
               'Expected %s.'], ...
               material, expected_dim, length(PT), shape_str);
    end
    for k = 1:expected_dim
        v = PT{k};
        if ~isnumeric(v) || ~isreal(v) || isempty(v) || ~isvector(v)
            error('SeaFreeze:badInput', ...
                  'PT{%d} must be a non-empty real numeric vector.', k);
        end
        if any(~isfinite(v(:)))
            error('SeaFreeze:badInput', ...
                  'PT{%d} contains non-finite values (Inf/NaN).', k);
        end
    end
    P = PT{1}; T = PT{2};
elseif isnumeric(PT)
    if ~ismatrix(PT) || size(PT,2) ~= expected_dim || isempty(PT)
        error('SeaFreeze:badInput', ...
              ['Scatter PT for material ''%s'' must be a non-empty N-by-%d ' ...
               'matrix; got size %s. Expected %s.'], ...
               material, expected_dim, mat2str(size(PT)), shape_str);
    end
    if ~isreal(PT) || any(~isfinite(PT(:)))
        error('SeaFreeze:badInput', ...
              'Scatter PT contains non-real or non-finite values (Inf/NaN).');
    end
    P = PT(:,1); T = PT(:,2);
else
    error('SeaFreeze:badInput', ...
          'PT must be a cell array or numeric matrix. Expected %s.', shape_str);
end

% Soft physical-range checks (warn rather than error: callers may legitimately
% probe extrapolation behaviour, and the spline returns NaN out of range).
if any(P(:) < 0)
    warning('SeaFreeze:negativePressure', ...
            'PT contains negative pressures; results will be NaN out of range.');
end
if any(T(:) <= 0)
    warning('SeaFreeze:nonPositiveT', ...
            'PT contains non-positive temperatures (T must be in Kelvin).');
end
if is_solution
    if iscell(PT), m = PT{3}; else, m = PT(:,3); end
    if any(m(:) < 0)
        error('SeaFreeze:badInput', 'Molality must be non-negative.');
    end
end
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end
