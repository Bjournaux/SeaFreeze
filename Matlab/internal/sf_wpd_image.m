function img = sf_wpd_image(idx, cols, cmiss, ctwo)
% SF_WPD_IMAGE  RGB image of a phase-index map (0 = not modelled,
% 1..n = phase, n+1 = two-phase region when ctwo is given).
img = zeros([size(idx) 3]);
for c = 1:3
    ch = cmiss(c) * ones(size(idx));
    for k = 1:size(cols, 1)
        ch(idx == k) = cols(k, c);
    end
    if nargin > 3
        ch(idx == size(cols, 1) + 1) = ctwo(c);
    end
    img(:, :, c) = ch;
end
end
