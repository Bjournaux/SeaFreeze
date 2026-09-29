function segs = sf_contour_segments(x, y, Z)
% SF_CONTOUR_SEGMENTS  Zero-level contour of Z(x,y) as a cell of [x y] segments.
%   Z is numel(y)-by-numel(x) (contourc convention).  NaN cells are excluded.
segs = {};
if all(isnan(Z(:))), return; end
C = contourc(x, y, Z, [0 0]);
i = 1;
while i <= size(C, 2)
    npt = C(2, i);
    if npt > 2, segs{end+1} = C(:, i+1:i+npt).'; end %#ok<AGROW>
    i = i + npt + 1;
end
end
