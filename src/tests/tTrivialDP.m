classdef tTrivialDP < matlab.unittest.TestCase
%tTrivialDP  Analytical unit tests for the core DP algorithm.
%
%   Uses two trivial problems with known optimal solutions:
%
%   Problem A (1 stage, direct-set dynamics):
%     x_next = u,  cost = u^2
%     StateFinal = {[0.9, 1.0]},  VFPenalty = 'rift'
%     Optimum: u* = 0.9,  totalCost = 0.81
%
%   Problem B (2 stages, additive dynamics):
%     x_next = x + u,  cost = u^2
%     StateFinal = {[1.0, 1.0]},  VFPenalty = 'rift'
%     Optimum: u*_seq = [0.5; 0.5],  totalCost = 0.5

    properties
        ProbA
        ProbB
    end

    methods (TestClassSetup)
        function setupProblems(tc)
            % 1 stage: x_next = u, cost = u^2, target [0.9, 1.0]
            tc.ProbA = DynaProg({0:0.1:1}, {0}, {[0.9, 1.0]}, {0:0.1:1}, 1, ...
                @trivial_sys, 'VFPenalty', 'rift', 'Display', 'off');
            tc.ProbA = run(tc.ProbA);

            % 2 stages: x_next = x+u, cost = u^2, target {1.0}
            tc.ProbB = DynaProg({0:0.05:2}, {0}, {[1.0, 1.0]}, {0:0.05:1}, 2, ...
                @trivial_additive_sys, 'VFPenalty', 'rift', 'Display', 'off');
            tc.ProbB = run(tc.ProbB);
        end
    end

    methods (Test)
        % --- Problem A: core results ---
        function testOptimalCost_SingleStage(tc)
            % Check the cost is equal to known value
            tc.verifyEqual(tc.ProbA.totalCost, 0.81, 'AbsTol', 1e-6)
        end

        function testOptimalControl_SingleStage(tc)
            % Check the CV profile is equal to known value
            tc.verifyEqual(tc.ProbA.ControlProfile{1}(1), 0.9, 'AbsTol', 1e-6)
        end

        function testFinalState_SingleStage(tc)
            % Check the final SV meets the target
            xf = tc.ProbA.StateProfile{1}(end);
            tc.verifyGreaterThanOrEqual(xf, 0.9)
            tc.verifyLessThanOrEqual(xf, 1.0)
        end

        % --- Problem B: test two-stage Bellman recursion ---
        function testOptimalCost_TwoStage(tc)
            % Check the cost is equal to known value
            tc.verifyEqual(tc.ProbB.totalCost, 0.5, 'AbsTol', 1e-6)
        end

        function testFinalState_TwoStage(tc)
            % Check the final state is at lower bound
            xf = tc.ProbB.StateProfile{1}(end);
            tc.verifyEqual(xf, 1.0, 'AbsTol', 0.001)
        end

        % --- Structural invariants (tested on Problem B) ---
        function testCostProfileSumMatchesTotalCost(tc)
            % Check that totalCost = sum of stage costs
            tc.verifyEqual(sum(tc.ProbB.CostProfile), tc.ProbB.totalCost, 'AbsTol', 1e-12)
        end

        function testCostProfileLength(tc)
            % Check that CostProfile has Nstages stage costs + 1 terminal cost
            tc.verifyNumElements(tc.ProbB.CostProfile, tc.ProbB.Nstages + 1)
        end

        function testStateProfileLength(tc)
            % Check that StateProfile has initial state + Nstages transitions
            tc.verifyNumElements(tc.ProbB.StateProfile{1}, tc.ProbB.Nstages + 1)
        end

        function testControlProfileLength(tc)
            % Check that ControlProfile has Nstages controls
            tc.verifyNumElements(tc.ProbB.ControlProfile{1}, tc.ProbB.Nstages)
        end
    end
end
