classdef test_SeaFreeze < matlab.unittest.TestCase
% Unit tests for SeaFreeze.m
% Mirrors the Python test suite at Python/seafreeze/test/test_seafreeze.py
%
% Run all tests:
%   results = runtests('test_SeaFreeze');
%   table(results)
%
% Notes:
%   - Tests must be run from the Matlab/ directory (or SeaFreeze_Gibbs.mat
%     must be on the MATLAB path) so that load() inside SeaFreeze.m succeeds.
%   - NaCl aqueous solution tests are omitted: the MATLAB SeaFreeze.m does
%     not expose a 'NaClaq' material case (unlike the Python getProp()).
%   - The shear modulus, Vp and Vs helper-function tests (Python
%     test_get_shear_mod_GPa_*, test_get_Vp_*, test_get_Vs_*) are covered
%     here by (a) verifying the formula directly against known values and
%     (b) verifying internal consistency of the SeaFreeze output, because
%     MATLAB computes these quantities inline rather than via separate
%     callable functions.

    % ------------------------------------------------------------------ %
    %  Shared expected values taken from the Python test suite
    % ------------------------------------------------------------------ %
    properties (Constant)
        % Ice VI at (900 MPa, 255 K)
        VI_P   = 900;
        VI_T   = 255;
        VI_rho = 1.356072490993616e+03;
        VI_Ks  = 1.832349756691741e+04;
        VI_shear = 7.303268592388283e+03;   % MPa
        VI_Vp  = 4.548954381485812e+03;     % m/s
        VI_Vs  = 2.320690281146717e+03;     % m/s

        % Ice III at (227 MPa, 244 K)
        III_P     = 227;
        III_T     = 244;
        III_shear = 3.9989e+03;   % MPa  (places=1)
        III_Vp    = 3.6242e+03;   % m/s  (places=1)
        III_Vs    = 1.8579e+03;   % m/s  (places=1)

        % Ice VII/X at (2500 MPa, 244 K)
        VII_P = 2500;
        VII_T = 244;
        VII_G = 1.899810e+06;   % J/kg  (places=0)

        % Shear-modulus parameters for each ice phase (same as SeaFreeze.m)
        SHEAR_Ih  = [3.1  -0.00462  0  -0.00657  1000  273.15];
        SHEAR_II  = [4.1   0.0175   0  -0.014    1100  273   ];
        SHEAR_III = [2.57  0.0175   0  -0.014    1100  273   ];
        SHEAR_V   = [2.57  0.0175   0  -0.014    1100  273   ];
        SHEAR_VI  = [2.57  0.0175   0  -0.014    1100  273   ];
        SHEAR_VII = [10    0.0033   0.000048 -0.014 1300  273  ];
    end

    % ================================================================== %
    %  Helper: compute shear modulus (GPa) from parameters, rho, T
    %  Mirrors Python seafreeze._get_shear_mod_GPa()
    % ================================================================== %
    methods (Static)
        function sm = shear_mod_GPa(parms, rho, T)
            sm = parms(1) ...
               + parms(2) .* (rho - parms(5)) ...
               + parms(3) .* (rho - parms(5)).^2 ...
               + parms(4) .* (T   - parms(6));
        end

        function Vp = get_Vp(shear_GPa, rho, Ks)
            % Mirrors Python seafreeze._get_Vp()
            % Ks in MPa, rho in kg/m^3, shear in GPa → Vp in m/s
            Vp = 1e3 * sqrt((Ks/1e3 + 4/3 .* shear_GPa) ./ rho / 1e-3);
        end

        function Vs = get_Vs(shear_GPa, rho)
            % Mirrors Python seafreeze._get_Vs()
            Vs = 1e3 * sqrt(shear_GPa ./ rho / 1e-3);
        end
    end

    % ================================================================== %
    %  _get_shear_mod_GPa equivalents  (Python: test_get_shear_mod_GPa_*)
    % ================================================================== %
    methods (Test)

        function test_shear_mod_singlept(testCase)
            % Mirrors test_get_shear_mod_GPa_singlept
            rho = testCase.VI_rho;
            T   = testCase.VI_T;
            sm  = testCase.shear_mod_GPa(testCase.SHEAR_VI, rho, T);
            testCase.verifyEqual(sm, 7.303268592388283, 'RelTol', 1e-10);
        end

        function test_shear_mod_multipt(testCase)
            % Mirrors test_get_shear_mod_GPa_multipt
            rho = [1356.072490993616; 1354.106807053618; 1352.057926458423];
            T   = [255; 265; 275];
            sm  = testCase.shear_mod_GPa(testCase.SHEAR_VI, rho, T);
            expected = [7.303268592388283; 7.128869123438319; 6.953013713022407];
            testCase.verifyEqual(sm, expected, 'RelTol', 1e-10);
        end

        function test_shear_mod_grid(testCase)
            % Mirrors test_get_shear_mod_GPa_grid
            % rho and T are already meshed (3×3 matrices)
            rho = [1356.072490993616  1353.307249697806  1350.440903858578;
                   1356.862314715232  1354.106807053618  1351.250712225762;
                   1357.649726085489  1354.903864412315  1352.057926458423];
            T   = repmat([255 265 275], 3, 1);   % T varies across columns
            sm  = testCase.shear_mod_GPa(testCase.SHEAR_VI, rho, T);
            expected = [7.303268592388283  7.114876869711612  6.924715817525114;
                        7.317090507516553  7.128869123438319  6.938887463950844;
                        7.330870206496061  7.142817627215505  6.953013713022400];
            testCase.verifyEqual(sm, expected, 'RelTol', 1e-10);
        end

        % ================================================================== %
        %  _get_Vp equivalents  (Python: test_get_Vp_*)
        % ================================================================== %

        function test_Vp_singlept(testCase)
            % Mirrors test_get_Vp_singlept
            smg = 7.303268592388283;
            rho = 1356.072490993616;
            Ks  = 18323.49756691741;
            Vp  = testCase.get_Vp(smg, rho, Ks);
            testCase.verifyEqual(Vp, 4548.954381485812, 'RelTol', 1e-10);
        end

        function test_Vp_multipt(testCase)
            % Mirrors test_get_Vp_multipt
            smg = [7.303268592388283; 7.128869123438319; 6.953013713022407];
            rho = [1356.072490993616; 1354.106807053618; 1352.057926458423];
            Ks  = [18323.49756691741; 18221.54159041028; 18115.56318084243];
            Vp  = testCase.get_Vp(smg, rho, Ks);
            expected = [4548.954381485812; 4525.042206621253; 4500.581390585054];
            testCase.verifyEqual(Vp, expected, 'RelTol', 1e-10);
        end

        function test_Vp_grid(testCase)
            % Mirrors test_get_Vp_grid
            smg = [7.303268592388283  7.114876869711612  6.924715817525114;
                   7.317090507516553  7.128869123438319  6.938887463950844;
                   7.330870206496061  7.142817627215505  6.953013713022400];
            rho = [1356.072490993616  1353.307249697806  1350.440903858578;
                   1356.862314715232  1354.106807053618  1351.250712225762;
                   1357.649726085489  1354.903864412315  1352.057926458423];
            Ks  = [18323.49756691741  18159.55544116658  17990.81337195430;
                   18385.04086077863  18221.54159041028  18053.25505449579;
                   18446.46522308364  18283.40144175383  18115.56318084243];
            Vp  = testCase.get_Vp(smg, rho, Ks);
            expected = [4548.954381485813  4519.791517059500  4489.896438587049;
                        4554.105836046733  4525.042206621253  4495.251116541151;
                        4559.235383167323  4530.269763950370  4500.581390585055];
            testCase.verifyEqual(Vp, expected, 'RelTol', 1e-10);
        end

        % ================================================================== %
        %  _get_Vs equivalents  (Python: test_get_Vs_*)
        % ================================================================== %

        function test_Vs_singlept(testCase)
            % Mirrors test_get_Vs_singlept
            smg = 7.303268592388283;
            rho = 1356.072490993616;
            Vs  = testCase.get_Vs(smg, rho);
            testCase.verifyEqual(Vs, 2320.690281146717, 'RelTol', 1e-10);
        end

        function test_Vs_multipt(testCase)
            % Mirrors test_get_Vs_multipt
            smg = [7.303268592388283; 7.128869123438319; 6.953013713022407];
            rho = [1356.072490993616; 1354.106807053618; 1352.057926458423];
            Vs  = testCase.get_Vs(smg, rho);
            expected = [2320.690281146717; 2294.477800749785; 2267.717197934159];
            testCase.verifyEqual(Vs, expected, 'RelTol', 1e-10);
        end

        function test_Vs_grid(testCase)
            % Mirrors test_get_Vs_grid
            smg = [7.303268592388283  7.114876869711612  6.924715817525114;
                   7.317090507516553  7.128869123438319  6.938887463950844;
                   7.330870206496061  7.142817627215505  6.953013713022400];
            rho = [1356.072490993616  1353.307249697806  1350.440903858578;
                   1356.862314715232  1354.106807053618  1351.250712225762;
                   1357.649726085489  1354.903864412315  1352.057926458423];
            Vs  = testCase.get_Vs(smg, rho);
            expected = [2320.690281146717  2292.901984107111  2264.452345731802;
                        2322.209103249143  2294.477800749785  2266.088955400825;
                        2323.720540877948  2296.045764570984  2267.717197934159];
            testCase.verifyEqual(Vs, expected, 'RelTol', 1e-10);
        end

        % ================================================================== %
        %  SF_getprop() - scatter single point  (Python: test_getProp_singlept)
        % ================================================================== %

        function test_getProp_singlept_VI(testCase)
            % Mirrors test_getProp_singlept
            % Scatter input: [P T] row vector
            out = SF_getprop([testCase.VI_P, testCase.VI_T], 'VI');
            testCase.verifyEqual(out.rho,   testCase.VI_rho,   'RelTol', 1e-5);
            testCase.verifyEqual(out.Ks,    testCase.VI_Ks,    'AbsTol', 1);
            testCase.verifyEqual(out.shear, testCase.VI_shear, 'RelTol', 1e-5);
            testCase.verifyEqual(out.Vp,    testCase.VI_Vp,    'AbsTol', 1);
            testCase.verifyEqual(out.Vs,    testCase.VI_Vs,    'RelTol', 1e-5);
        end

        function test_getProp_singlept_III_shear(testCase)
            % Mirrors test_getProp_shear_mod_parms_point
            out = SF_getprop([testCase.III_P, testCase.III_T], 'III');
            testCase.verifyEqual(out.shear, testCase.III_shear, 'AbsTol', 10);
            testCase.verifyEqual(out.Vp,    testCase.III_Vp,    'AbsTol', 10);
            testCase.verifyEqual(out.Vs,    testCase.III_Vs,    'AbsTol', 10);
        end

        function test_getProp_VII_G(testCase)
            % Mirrors test_getProp_VII_pt — Gibbs energy of ice VII/X
            out = SF_getprop([testCase.VII_P, testCase.VII_T], 'VII_X_French');
            testCase.verifyEqual(out.G, testCase.VII_G, 'AbsTol', 1);
        end

        % ================================================================== %
        %  SF_getprop() - scatter multiple points
        % ================================================================== %

        function test_getProp_scatter_VI_multipt(testCase)
            % Three ice VI points used throughout the shear / Vp / Vs tests
            PT = [900 255; 910 265; 920 275];
            out = SF_getprop(PT, 'VI');
            expected_rho = [1356.072490993616; 1354.106807053618; 1352.057926458423];
            testCase.verifyEqual(out.rho, expected_rho, 'RelTol', 1e-5);
        end

        function test_getProp_scatter_VI_shear_multipt(testCase)
            % Shear for three ice VI scatter points
            PT = [900 255; 910 265; 920 275];
            out = SF_getprop(PT, 'VI');
            expected_shear = [7.303268592388283; 7.128869123438319; 6.953013713022407] * 1e3;
            testCase.verifyEqual(out.shear, expected_shear, 'RelTol', 1e-5);
        end

        function test_getProp_scatter_VI_Vp_multipt(testCase)
            PT = [900 255; 910 265; 920 275];
            out = SF_getprop(PT, 'VI');
            expected_Vp = [4548.954381485812; 4525.042206621253; 4500.581390585054];
            testCase.verifyEqual(out.Vp, expected_Vp, 'RelTol', 1e-4);
        end

        function test_getProp_scatter_VI_Vs_multipt(testCase)
            PT = [900 255; 910 265; 920 275];
            out = SF_getprop(PT, 'VI');
            expected_Vs = [2320.690281146717; 2294.477800749785; 2267.717197934159];
            testCase.verifyEqual(out.Vs, expected_Vs, 'RelTol', 1e-5);
        end

        % ================================================================== %
        %  SF_getprop() - gridded input  (Python: test_getProp with grid PT)
        % ================================================================== %

        function test_getProp_grid_VI(testCase)
            % Mirrors the grid variant of test_get_shear_mod_GPa_grid
            % P: 900..920 (3 pts), T: 255..275 (3 pts) → 3×3 output
            P = (900:10:920)';
            T = (255:10:275)';
            out = SF_getprop({P, T}, 'VI');
            % Corner (P=900, T=255) must match the single-point result
            testCase.verifyEqual(out.rho(1,1),   testCase.VI_rho,   'RelTol', 1e-5);
            testCase.verifyEqual(out.shear(1,1), testCase.VI_shear, 'RelTol', 1e-5);
            testCase.verifyEqual(out.Vp(1,1),    testCase.VI_Vp,    'AbsTol', 1);
            testCase.verifyEqual(out.Vs(1,1),    testCase.VI_Vs,    'RelTol', 1e-5);
        end

        function test_getProp_grid_rho_matrix(testCase)
            % Full 3×3 rho grid for ice VI (P: 900..920, T: 255..275)
            P = (900:10:920)';
            T = (255:10:275)';
            out = SF_getprop({P, T}, 'VI');
            % Row 3, Col 3 → (920 MPa, 275 K)
            expected_rho_33 = 1352.057926458423;
            testCase.verifyEqual(out.rho(3,3), expected_rho_33, 'RelTol', 1e-5);
        end

        function test_getProp_grid_shear_matrix(testCase)
            % Full 3×3 shear grid for ice VI
            P = (900:10:920)';
            T = (255:10:275)';
            out = SF_getprop({P, T}, 'VI');
            expected_shear = [7.303268592388283  7.114876869711612  6.924715817525114;
                              7.317090507516553  7.128869123438319  6.938887463950844;
                              7.330870206496061  7.142817627215505  6.953013713022400] * 1e3;
            testCase.verifyEqual(out.shear, expected_shear, 'RelTol', 1e-5);
        end

        function test_getProp_grid_Vp_matrix(testCase)
            % Full 3×3 Vp grid for ice VI
            P = (900:10:920)';
            T = (255:10:275)';
            out = SF_getprop({P, T}, 'VI');
            expected_Vp = [4548.954381485813  4519.791517059500  4489.896438587049;
                           4554.105836046733  4525.042206621253  4495.251116541151;
                           4559.235383167323  4530.269763950370  4500.581390585055];
            testCase.verifyEqual(out.Vp, expected_Vp, 'RelTol', 1e-4);
        end

        function test_getProp_grid_Vs_matrix(testCase)
            % Full 3×3 Vs grid for ice VI
            P = (900:10:920)';
            T = (255:10:275)';
            out = SF_getprop({P, T}, 'VI');
            expected_Vs = [2320.690281146717  2292.901984107111  2264.452345731802;
                           2322.209103249143  2294.477800749785  2266.088955400825;
                           2323.720540877948  2296.045764570984  2267.717197934159];
            testCase.verifyEqual(out.Vs, expected_Vs, 'RelTol', 1e-5);
        end

        function test_getProp_grid_output_size(testCase)
            % Grid output shape must equal [numel(P), numel(T)]
            P = (400:2:500)';   % 51 values
            T = (220:0.5:250)'; % 61 values
            out = SF_getprop({P, T}, 'V');
            testCase.verifySize(out.rho, [numel(P), numel(T)]);
        end

        % ================================================================== %
        %  SF_getprop() - ice Ih  (scatter single point)
        % ================================================================== %

        function test_getProp_Ih_density_range(testCase)
            % Ice Ih near melting point - density slightly below 1000 kg/m^3
            out = SF_getprop([1, 265], 'Ih');
            testCase.verifyGreaterThan(out.rho, 800);
            testCase.verifyLessThan(out.rho, 950);
        end

        function test_getProp_Ih_has_shear(testCase)
            % Ice Ih is a solid - must have shear modulus and wave velocities
            out = SF_getprop([1, 265], 'Ih');
            testCase.verifyTrue(isfield(out, 'shear'));
            testCase.verifyGreaterThan(out.shear, 0);
            testCase.verifyGreaterThan(out.Vp, 0);
            testCase.verifyGreaterThan(out.Vs, 0);
        end

        % ================================================================== %
        %  SF_getprop() - liquid water phases (no shear)
        % ================================================================== %

        function test_getProp_water1_no_shear(testCase)
            % Liquid water: SeaFreeze.m does not set shear for water phases
            out = SF_getprop([500, 300], 'water1');
            testCase.verifyFalse(isfield(out, 'shear'));
        end

        function test_getProp_water1_positive_density(testCase)
            out = SF_getprop([500, 300], 'water1');
            testCase.verifyGreaterThan(out.rho, 500);
        end

        function test_getProp_water1_grid(testCase)
            P = (100:200:500)';
            T = (250:25:300)';
            out = SF_getprop({P, T}, 'water1');
            testCase.verifySize(out.rho, [numel(P), numel(T)]);
            testCase.verifyGreaterThan(min(out.rho(:)), 500);
        end

        % ================================================================== %
        %  Thermodynamic consistency checks
        % ================================================================== %

        function test_enthalpy_consistency(testCase)
            % H = G + T*S  (from definition of Gibbs energy)
            out = SF_getprop([testCase.VI_P, testCase.VI_T], 'VI');
            H_calc = out.G + testCase.VI_T * out.S;
            testCase.verifyEqual(out.H, H_calc, 'RelTol', 1e-8);
        end

        function test_internal_energy_consistency(testCase)
            % U = H - P/rho = G + T*S - P/rho
            % (P in MPa → multiply by 1e6 for Pa, rho in kg/m^3)
            out = SF_getprop([testCase.VI_P, testCase.VI_T], 'VI');
            U_calc = out.G + testCase.VI_T * out.S ...
                     - testCase.VI_P * 1e6 / out.rho;
            testCase.verifyEqual(out.U, U_calc, 'RelTol', 1e-8);
        end

        function test_Ks_from_bulk_sound_speed(testCase)
            % Ks = rho * vel^2 / 1e6  (vel = bulk sound speed in m/s)
            out = SF_getprop([testCase.VI_P, testCase.VI_T], 'VI');
            Ks_calc = out.rho * out.vel^2 / 1e6;
            testCase.verifyEqual(out.Ks, Ks_calc, 'RelTol', 1e-8);
        end

        function test_shear_formula_VI_singlept(testCase)
            % Verify SeaFreeze.m shear is consistent with the published formula
            % shear [GPa] = a + b*(rho-rho0) + c*(rho-rho0)^2 + d*(T-T0)
            out = SF_getprop([testCase.VI_P, testCase.VI_T], 'VI');
            sm_GPa = testCase.shear_mod_GPa(testCase.SHEAR_VI, out.rho, testCase.VI_T);
            testCase.verifyEqual(out.shear, sm_GPa * 1e3, 'RelTol', 1e-10);
        end

        function test_Vp_formula_VI_singlept(testCase)
            % Vp = 1e3 * sqrt((Ks/1e3 + 4/3*shear_GPa) / (rho*1e-3))
            out = SF_getprop([testCase.VI_P, testCase.VI_T], 'VI');
            shear_GPa = out.shear / 1e3;
            Vp_calc = 1e3 * sqrt((out.Ks/1e3 + 4/3 * shear_GPa) / (out.rho * 1e-3));
            testCase.verifyEqual(out.Vp, Vp_calc, 'RelTol', 1e-10);
        end

        function test_Vs_formula_VI_singlept(testCase)
            % Vs = 1e3 * sqrt(shear_GPa / (rho*1e-3))
            out = SF_getprop([testCase.VI_P, testCase.VI_T], 'VI');
            shear_GPa = out.shear / 1e3;
            Vs_calc = 1e3 * sqrt(shear_GPa / (out.rho * 1e-3));
            testCase.verifyEqual(out.Vs, Vs_calc, 'RelTol', 1e-10);
        end

        % ================================================================== %
        %  Scatter vs. grid consistency
        % ================================================================== %

        function test_scatter_and_grid_agree_VI(testCase)
            % A single-point scatter result must match the corresponding
            % element of a grid result.
            P0 = 900; T0 = 255;
            P  = (900:10:920)';
            T  = (255:10:275)';

            out_scatter = SF_getprop([P0, T0], 'VI');
            out_grid    = SF_getprop({P, T}, 'VI');

            testCase.verifyEqual(out_scatter.rho,   out_grid.rho(1,1),   'RelTol', 1e-8);
            testCase.verifyEqual(out_scatter.Ks,    out_grid.Ks(1,1),    'RelTol', 1e-8);
            testCase.verifyEqual(out_scatter.shear, out_grid.shear(1,1), 'RelTol', 1e-8);
            testCase.verifyEqual(out_scatter.Vp,    out_grid.Vp(1,1),    'RelTol', 1e-8);
            testCase.verifyEqual(out_scatter.Vs,    out_grid.Vs(1,1),    'RelTol', 1e-8);
        end

    end
end
