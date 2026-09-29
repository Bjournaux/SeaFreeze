function [cols, cmiss, ctwo] = sf_wpd_colors(n)
% SF_WPD_COLORS  Field colours for SF_WPD_PT / SF_WPD_rhoT (light tints of a
% fixed categorical order: fluid, Ih, II, III, V, VI, VII/X), plus the colour
% of 'not modelled' cells and of two-phase regions.
base = [0.165 0.471 0.839; 0.922 0.408 0.204; 0.106 0.686 0.478; 0.929 0.631 0.000; ...
        0.910 0.482 0.643; 0.000 0.514 0.000; 0.290 0.227 0.655];
a = 0.30;
cols = 1 - a * (1 - base(1:n, :));
cmiss = [0.93 0.93 0.92];
ctwo  = [0.86 0.86 0.84];
end
