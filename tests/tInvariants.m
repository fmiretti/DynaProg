classdef tInvariants < matlab.unittest.TestCase
%tInvariants  Mathematical and structural equivalence tests.
%
%   Tests invariants that must hold regardless of implementation flags:
%     - SafeMode=true gives the same result as SafeMode=false
%     - Split-model gives the same result as a single model
%     - An infeasible problem returns totalCost = Inf

    properties
        ProbSafeModeOff
        ProbSafeModeOn
    end

    methods (TestClassSetup)
        function setupSafeModeProblems(tc)
            tc.ProbSafeModeOff = tc.solveProbA('SafeMode', false);
            tc.ProbSafeModeOn  = tc.solveProbA('SafeMode', true);
        end
    end

    methods (Test)
        % --- SafeMode equivalence ---
        function testSafeMode_SameTotalCost(tc)
            tc.verifyEqual(tc.ProbSafeModeOff.totalCost, tc.ProbSafeModeOn.totalCost, 'AbsTol', 1e-10)
        end

        function testSafeMode_SameControlProfile(tc)
            tc.verifyEqual(tc.ProbSafeModeOff.ControlProfile{1}, tc.ProbSafeModeOn.ControlProfile{1}, 'AbsTol', 1e-10)
        end

        function testSafeMode_SameStateProfile(tc)
            tc.verifyEqual(tc.ProbSafeModeOff.StateProfile{1}, tc.ProbSafeModeOn.StateProfile{1}, 'AbsTol', 1e-10)
        end

        % --- Split-model equivalence ---
        function testSplitModel_SameTotalCostAsSingleModel(tc)
            % Both models compute x_next=u, cost=u^2; results must be identical.
            % ExogenousInput is a dummy zeros vector so model_wrapper calls
            % SysNameInt with 4 args (required for the split-model int signature).
            exo = zeros(1, 1);
            probSingle = DynaProg({0:0.1:1}, {0}, {[0.9, 1.0]}, {0:0.1:1}, 1, ...
                @trivial_sys, ...
                'ExogenousInput', exo, 'Display', 'off');
            probSingle = run(probSingle);

            probSplit = DynaProg({0:0.1:1}, {0}, {[0.9, 1.0]}, {0:0.1:1}, 1, ...
                @trivial_ext, @trivial_int, ...
                'ExogenousInput', exo, 'Display', 'off');
            probSplit = run(probSplit);

            tc.verifyEqual(probSplit.totalCost, probSingle.totalCost, 'AbsTol', 1e-10)
        end

        function testSplitModel_SameControlProfile(tc)
            exo = zeros(1, 1);
            probSingle = DynaProg({0:0.1:1}, {0}, {[0.9, 1.0]}, {0:0.1:1}, 1, ...
                @trivial_sys, ...
                'ExogenousInput', exo, 'Display', 'off');
            probSingle = run(probSingle);

            probSplit = DynaProg({0:0.1:1}, {0}, {[0.9, 1.0]}, {0:0.1:1}, 1, ...
                @trivial_ext, @trivial_int, ...
                'ExogenousInput', exo, 'Display', 'off');
            probSplit = run(probSplit);

            tc.verifyEqual(probSplit.ControlProfile{1}, probSingle.ControlProfile{1}, 'AbsTol', 1e-10)
        end

        % --- Infeasible problem ---
        function testInfeasibleProblem_ReturnsTotalCostInf(tc)
            % Target xf \in [2.0, 3.0] is outside ControlGrid [0, 1]: no control can reach it.
            % Display='warn' keeps DisplayWarnings=true so failedBackward is detected;
            prob = DynaProg({0:0.1:1}, {0}, {[2.0, 3.0]}, {0:0.1:1}, 1, ...
                @trivial_sys, 'Display', 'warn');
            % Suppress MATLAB warnings to avoid clutter
            ws = warning('off', 'all');
            prob = run(prob);
            warning(ws);
            tc.verifyEqual(prob.totalCost, Inf)
        end
    end

    methods (Access = private)
        function prob = solveProbA(~, varargin)
            % Problem A with optional extra NVPs specified in varargin
            prob = DynaProg({0:0.1:1}, {0}, {[0.9, 1.0]}, {0:0.1:1}, 1, ...
                @trivial_sys, 'Display', 'off', varargin{:});
            prob = run(prob);
        end
    end
end
