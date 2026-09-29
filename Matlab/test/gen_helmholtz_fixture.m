% gen_helmholtz_fixture.m — build a test Helmholtz spline F(rho,T) for water.
%
% Samples the analytic IAPWS-95 Helmholtz energy (legacy/LocalBasisFunction/
% IAPWS95.m, rho-T mode) on a liquid-region grid and fits a tensor B-spline
% with spdft.  The result is a stand-in for a production F(rho,T) spline so
% that fnFval / getProp's Helmholtz path can be tested end to end.
%
% Output: test/fixtures/water_F_test.mat (v7, variable 'sp'), containing
%   sp.form/knots/number/order/dim/coefs   B-form, knots{1}=rho, knots{2}=T
%   sp.eos    = 'F_rhoT'                   Helmholtz-spline flag
%   sp.MW     = 0.018015268                kg/mol
%   sp.source = description string
%
% Run from Matlab/ (MATLAB or Octave):  run('test/gen_helmholtz_fixture.m')

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root, 'internal'));
addpath(fullfile(root, 'internal', 'compat'));
addpath(fullfile(root, 'legacy'));
addpath(fullfile(root, 'legacy', 'LocalBasisFunction'));

% spdft stamps sp.revision with datetime, which base Octave lacks.
if exist('OCTAVE_VERSION', 'builtin') && ~exist('datetime')
    shimdir = tempname(); mkdir(shimdir);
    fid = fopen(fullfile(shimdir, 'datetime.m'), 'w');
    fprintf(fid, 'function d = datetime()\nd = datestr(now);\nend\n');
    fclose(fid);
    addpath(shimdir);
end

rho = 900:4:1400;     % kg/m^3
T   = 240:2:500;      % K

fprintf('Sampling IAPWS-95 on %d x %d (rho,T) grid...\n', numel(rho), numel(T));
r = IAPWS95({rho, T}, 'rho');
A = real(r.A);        % J/kg
P = real(r.P);        % MPa

% IAPWS-95 is unphysical beyond the liquid spinodal (P blows up to 1e6 MPa
% at low rho / low T).  Fit only the mechanically stable region; the
% regularization fills the masked corner smoothly.
dPdr = [diff(P,1,1); diff(P(end-1:end,:),1,1)];
ok   = P > -60 & dPdr > 0 & real(r.Kt) > 0;
mask = ones(size(A)); mask(~ok) = NaN;

opts.Xc     = {linspace(rho(1), rho(end), 30), linspace(T(1), T(end), 30)};
opts.ordr   = [6 6];
opts.mdrv   = [4 4];
opts.lam    = [1e-4 1e-4];
opts.RegFac = [2 2];
opts.mask   = mask;

fprintf('Fitting spline...\n');
sp = spdft({rho, T}, A, 1, opts);
sp = rmfield(sp, 'Data');     % keep the fixture small
sp.eos    = 'F_rhoT';
sp.MW     = 0.018015268;
sp.source = 'IAPWS-95 A(rho,T) fit by gen_helmholtz_fixture.m (test fixture only)';

Afit = sp_val(sp, {rho, T});
fprintf('max |A_fit - A_IAPWS| = %.3g J/kg (fitted region)\n', max(abs(Afit(ok) - A(ok))));

outdir = fullfile(here, 'fixtures');
if ~exist(outdir, 'dir'), mkdir(outdir); end
save(fullfile(outdir, 'water_F_test.mat'), 'sp', '-v7');
fprintf('Wrote %s\n', fullfile(outdir, 'water_F_test.mat'));
