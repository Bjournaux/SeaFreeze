% water3_tutorial.m — SeaFreeze tutorial: the Helmholtz fluid 'water3',
% density input, and liquid-vapour-ice coexistence.  Runs under MATLAB and
% GNU Octave.  Every example of the "water3" section of Matlab/README.md,
% plus a four-panel figure saved as water3_tutorial_matlab.png next to this
% file.
%
%   run('examples/water3_tutorial.m')      % from the Matlab/ folder

here = fileparts(mfilename('fullpath'));
addpath(genpath(fileparts(here)));

%% 1. water3 at (P,T): same call and outputs as every other SeaFreeze phase
out = SF_getprop([0.101325 298.15; 100 300; 1000 350], 'water3');
fprintf('1. water3: rho = %s kg/m3, Cp = %s J/kg/K, vel = %s m/s\n', ...
        mat2str(out.rho', 7), mat2str(out.Cp', 5), mat2str(out.vel', 5));
g3 = SF_getprop({[0.1 50 500], [280 300 320]}, 'water3', {'rho','alpha'});   % grid: P down, T across
g1 = SF_getprop({[0.1 50 500], [280 300 320]}, 'water1', {'rho','alpha'});
fprintf('   grid rho water3 - water1 (kg/m3):\n'); disp(round(g3.rho - g1.rho, 3));

%% 2. One fluid, two branches: the stable phase is the lower Gibbs energy
T = linspace(280, 500, 221)';
stable = SF_getprop([repmat(0.101325, numel(T), 1) T], 'water3', {'rho','G'});
liquid = SF_getprop([repmat(0.101325, numel(T), 1) T], 'water3', {'rho','G'}, 'branch', 'liquid');
Tb = T(find(stable.rho < 100, 1));
fprintf('2. at 0.101325 MPa the stable branch switches to vapour at %.0f K (T_sat 373.124 K)\n', Tb);

%% 3. Density-temperature input, for water3 AND the Gibbs splines
w3 = SF_getprop([1000 300; 1100 300], 'water3', {'P','Cp'}, 'input', 'rhoT');   % direct
w1 = SF_getprop([1000 300; 1100 300], 'water1', {'P','Cp'}, 'input', 'rhoT');   % via SF_rho2P
fprintf('3. P at 1000/1100 kg/m3, 300 K: water3 %s MPa, water1 %s MPa\n', mat2str(w3.P', 6), mat2str(w1.P', 6));
ice   = SF_getprop([1330 260], 'VI', {'P','Vp','Vs'}, 'input', 'rhoT');
brine = SF_getprop([1050 300 1.0], 'NaClaq', {'P','muw'}, 'input', 'rhoT');
fprintf('   ice VI at 1330 kg/m3, 260 K: P = %.1f MPa, Vp = %.0f m/s; NaCl(aq) 1 mol/kg, 1050 kg/m3, 300 K: P = %.2f MPa\n', ...
        ice.P, ice.Vp, brine.P);
fprintf('   SF_rho2P (water3): %s MPa\n', mat2str(SF_rho2P([997.047 1100], [298.15 300], 'water3')', 6));

%% 4. Phase equilibria with water3 as the liquid
ph = SF_WhichPhase({[0.1 300 800], [260 270 276 300]}, 'liquid', 'water3');
fprintf('4. SF_WhichPhase (0 = liquid, 1 = Ih, 6 = VI):\n'); disp(ph);
% Ih-liquid line down to the triple-point pressure (default grid starts at 0.1 MPa)
line = SF_PhaseLines('Ih', 'water3', 'P', logspace(log10(6.2e-4), log10(209), 300), 'T', 250:0.02:273.3);

%% 5. Saturation and sublimation curves
sat  = SF_coexistence('saturation', linspace(230, 646.5, 200));   % < 273.16 K: supercooled liquid
Tsub = linspace(170, 273.16, 150);
sub  = SF_coexistence('sublimation', Tsub);                        % < 230 K: dilute-vapour extension (warns once)
s1 = SF_coexistence('saturation', 373.124); s2 = SF_coexistence('sublimation', 250);
fprintf('5. p_sat(373.124 K) = %.6f MPa; p_sub(250 K) = %.3f Pa\n', s1.P, s2.P * 1e6);

%% Figure
fig = figure('Visible', 'off', 'Position', [100 100 1200 900]);
c3 = [0.165 0.471 0.839]; c1 = [0.922 0.408 0.204]; cv = [0.106 0.686 0.478];

subplot(2, 2, 1);
semilogy(T, stable.rho, '-', 'Color', c3, 'LineWidth', 2); hold on
semilogy(T, liquid.rho, '--', 'Color', c1, 'LineWidth', 1.5);
xlabel('Temperature (K)'); ylabel('Density (kg/m^3)'); grid on
title('Boiling at 0.101325 MPa: one EOS, two branches');
legend({'''stable'' (default)', '''liquid'' (metastable above T_b)'}, 'Location', 'east');

subplot(2, 2, 2);
rho = linspace(950, 1250, 61);
p3a = SF_getprop({rho, 280}, 'water3', 'P', 'input', 'rhoT');
p1a = SF_getprop({rho, 280}, 'water1', 'P', 'input', 'rhoT');
p3b = SF_getprop({rho, 350}, 'water3', 'P', 'input', 'rhoT');
p1b = SF_getprop({rho, 350}, 'water1', 'P', 'input', 'rhoT');
plot(rho, p3a.P, '-', 'Color', c3, 'LineWidth', 2); hold on
plot(rho, p1a.P, '-', 'Color', c1, 'LineWidth', 1.2);
plot(rho, p3b.P, '--', 'Color', c3, 'LineWidth', 2);
plot(rho, p1b.P, '--', 'Color', c1, 'LineWidth', 1.2);
xlabel('Density (kg/m^3)'); ylabel('Pressure (MPa)'); grid on
title('(\rho,T) input: isotherms P(\rho)');
legend({'water3, 280 K', 'water1 (via SF\_rho2P), 280 K', 'water3, 350 K', 'water1, 350 K'}, 'Location', 'northwest');

subplot(2, 2, 3);
semilogy(sat.T, sat.P * 1e6, '-', 'Color', c3, 'LineWidth', 2); hold on
in = Tsub >= 230;
semilogy(Tsub(in), sub.P(in) * 1e6, '-', 'Color', cv, 'LineWidth', 2);
semilogy(Tsub(~in), sub.P(~in) * 1e6, '--', 'Color', cv, 'LineWidth', 1.5);
[Ps, io] = sort(line.P);
semilogy(line.T(io), Ps * 1e6, '-', 'Color', c1, 'LineWidth', 2);
semilogy(273.16, 611.657, 'ks', 'MarkerFaceColor', 'k');
ylim([1e-3 3e8]); xlabel('Temperature (K)'); ylabel('Pressure (Pa)'); grid on
title('Triple point: melting, boiling, sublimation');
legend({'liquid-vapour', 'ice Ih-vapour', 'dilute-vapour extension', 'ice Ih-liquid'}, 'Location', 'southeast');

ax = subplot(2, 2, 4);
SF_WPD('ax', ax, 'liquid', 'water3', 'meta', 'none');
title(ax, 'SF\_WPD(''liquid'',''water3'')');

print(fig, fullfile(here, 'water3_tutorial_matlab.png'), '-dpng', '-r130');
fprintf('saved %s\n', fullfile(here, 'water3_tutorial_matlab.png'));
