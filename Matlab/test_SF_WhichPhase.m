classdef test_SF_WhichPhase < matlab.unittest.TestCase
% Unit tests for SF_WhichPhase.m
% Mirrors the Python test suite at Python/seafreeze/test/test_whichphase.py
%
% Run all tests:
%   results = runtests('test_SF_WhichPhase');
%   table(results)
%
% Requirements:
%   - MATLAB Curve Fitting Toolbox (SF_WhichPhase.m uses fnval / fnder).
%     Tests are skipped automatically when the toolbox is absent.
%   - SeaFreeze_Gibbs.mat must be on the MATLAB path (run from Matlab/).
%
% Phase code mapping (same as Python):
%   0 = liquid water    3 = ice III
%   1 = ice Ih          5 = ice V
%   2 = ice II          6 = ice VI
%   NaN = outside all parametrisation ranges
%
% Omitted tests (no MATLAB equivalent):
%   test_NaCl_grid, test_NaCl_singlePt, test_NaCl_singlePt2
%   — SF_WhichPhase.m does not include NaCl aqueous solutions.

    % ------------------------------------------------------------------ %
    %  Skip all tests when the Curve Fitting Toolbox is unavailable
    % ------------------------------------------------------------------ %
    methods (TestClassSetup)
        function checkToolbox(testCase)
            if ~license('test', 'curve_fitting_toolbox')
                testCase.assumeTrue(false, ...
                    'Curve Fitting Toolbox not available; skipping SF_WhichPhase tests.');
            end
        end
    end

    % ================================================================== %
    %  Scatter – single point  (Python: test_scatter_singlePt)
    % ================================================================== %
    methods (Test)

        function test_scatter_singlePt_iceVI(testCase)
            % (2000 MPa, 334 K) → ice VI
            out = SF_WhichPhase([2000, 334]);
            testCase.verifyEqual(out, 6);
        end

        function test_scatter_singlePt_water(testCase)
            % (500 MPa, 300 K) → liquid water
            out = SF_WhichPhase([500, 300]);
            testCase.verifyEqual(out, 0);
        end

        % ================================================================== %
        %  Scatter – multiple points  (Python: test_scatter_multiPt)
        % ================================================================== %

        function test_scatter_multiPt(testCase)
            % Mirrors test_scatter_multiPt
            % (100,200)→Ih, (400,250)→V, (1000,300)→VI
            % SF_WhichPhase returns a [1×n] row when given an [n×2] scatter matrix.
            PT = [100 200; 400 250; 1000 300];
            out = SF_WhichPhase(PT);
            testCase.verifyEqual(out(:)', [1 5 6]);
        end

        % ================================================================== %
        %  Scatter – all phases  (Python: test_phases)
        % ================================================================== %

        function test_phases_water(testCase)
            out = SF_WhichPhase([500, 300]);
            testCase.verifyEqual(out, 0);
        end

        function test_phases_Ih(testCase)
            out = SF_WhichPhase([1, 200]);
            testCase.verifyEqual(out, 1);
        end

        function test_phases_II(testCase)
            out = SF_WhichPhase([300, 200]);
            testCase.verifyEqual(out, 2);
        end

        function test_phases_III(testCase)
            out = SF_WhichPhase([300, 250]);
            testCase.verifyEqual(out, 3);
        end

        function test_phases_V(testCase)
            out = SF_WhichPhase([500, 250]);
            testCase.verifyEqual(out, 5);
        end

        function test_phases_VI(testCase)
            out = SF_WhichPhase([800, 250]);
            testCase.verifyEqual(out, 6);
        end

        function test_phases_VI_extremeP(testCase)
            % Extreme pressure at low T → still ice VI (Python: PT[6]=(2000,334))
            out = SF_WhichPhase([2000, 334]);
            testCase.verifyEqual(out, 6);
        end

        % ================================================================== %
        %  Scatter – extrapolation  (Python: test_scatter_all_extrapolations)
        % ================================================================== %

        function test_scatter_extrapolation_badP(testCase)
            % P=1e6 MPa is far outside all parametrisations → NaN
            out = SF_WhichPhase([1e6, 250]);
            testCase.verifyTrue(isnan(out), ...
                'Expected NaN for out-of-range pressure');
        end

        function test_scatter_extrapolation_badT(testCase)
            % T=5000 K is far outside all parametrisations → NaN
            out = SF_WhichPhase([400, 5000]);
            testCase.verifyTrue(isnan(out), ...
                'Expected NaN for out-of-range temperature');
        end

        function test_scatter_extrapolation_badPT(testCase)
            % Both P and T out of range → NaN
            out = SF_WhichPhase([2000, 2000]);
            testCase.verifyTrue(isnan(out), ...
                'Expected NaN for out-of-range P and T');
        end

        % ================================================================== %
        %  Grid  (Python: test_grid)
        % ================================================================== %

        function test_grid_output_size(testCase)
            % Grid output must have shape [numel(P), numel(T)]
            P = 0:200:1000;   % 6 values
            T = 200:50:350;   % 4 values
            out = SF_WhichPhase({P, T});
            testCase.verifySize(out, [numel(P), numel(T)]);
        end

        function test_grid_corner_Ih(testCase)
            % (P=0, T=200) → ice Ih
            P = 0:200:1000;
            T = 200:50:350;
            out = SF_WhichPhase({P, T});
            testCase.verifyEqual(out(1,1), 1);
        end

        function test_grid_spot_checks(testCase)
            % Mirrors Python test_grid expected array:
            %   P=[0,200,400,600,800,1000], T=[200,250,300,350]
            %   exp=[[1,1,0,0],[2,1,0,0],[2,5,0,0],[2,5,0,0],[6,6,0,0],[6,6,6,0]]
            % 0 in the Python array denotes entries where fnval returns 0
            % (all splines return 0 = out of range). MATLAB SF_WhichPhase
            % sets those to NaN instead, so we test only the non-zero entries.
            P = 0:200:1000;
            T = 200:50:350;
            out = SF_WhichPhase({P, T});

            % P=0,   T=200 → Ih (1)
            testCase.verifyEqual(out(1,1), 1);
            % P=0,   T=250 → Ih (1)
            testCase.verifyEqual(out(1,2), 1);
            % P=200, T=200 → II (2)
            testCase.verifyEqual(out(2,1), 2);
            % P=400, T=250 → V  (5)
            testCase.verifyEqual(out(3,2), 5);
            % P=800, T=200 → VI (6)
            testCase.verifyEqual(out(5,1), 6);
            % P=800, T=250 → VI (6)
            testCase.verifyEqual(out(5,2), 6);
            % P=1000,T=300 → VI (6)
            testCase.verifyEqual(out(6,3), 6);
        end

        % ================================================================== %
        %  Phase code validity
        % ================================================================== %

        function test_output_values_are_valid_phase_codes(testCase)
            % All returned values must belong to {0,1,2,3,5,6} or be NaN
            P = 0:100:1000;
            T = 220:10:350;
            out = SF_WhichPhase({P, T});
            valid = {0, 1, 2, 3, 5, 6};
            flat = out(:);
            for k = 1:numel(flat)
                v = flat(k);
                if ~isnan(v)
                    testCase.verifyTrue(any(cellfun(@(x) x==v, valid)), ...
                        sprintf('Unexpected phase code %g at index %d', v, k));
                end
            end
        end

    end
end
