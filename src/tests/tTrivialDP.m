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

    methods (Test)
        % --- Problem A: core results ---
        function testOptimalCost_SingleStage(tc)
            % Check the cost is equal to known value
            prob = tc.solveProbA();
            tc.verifyEqual(prob.totalCost, 0.81, 'AbsTol', 1e-6)
        end

        function testOptimalControl_SingleStage(tc)
            % Check the CV profile is equal to known value
            prob = tc.solveProbA();
            tc.verifyEqual(prob.ControlProfile{1}(1), 0.9, 'AbsTol', 1e-6)
        end

        function testFinalState_SingleStage(tc)
            % Check the final SV meets the target
            prob = tc.solveProbA();
            xf = prob.StateProfile{1}(end);
            tc.verifyGreaterThanOrEqual(xf, 0.9)
            tc.verifyLessThanOrEqual(xf, 1.0)
        end

        % --- Problem B: test two-stage Bellman recursion ---
        function testOptimalCost_TwoStage(tc)
            % Check the cost is equal to known value
            prob = tc.solveProbB();
            tc.verifyEqual(prob.totalCost, 0.5, 'AbsTol', 1e-6)
        end

        function testFinalState_TwoStage(tc)
            % Check the final state is at lower bound
            prob = tc.solveProbB();
            xf = prob.StateProfile{1}(end);
            tc.verifyEqual(xf, 1.0, 'AbsTol', 0.001)
        end

        % --- Structural invariants (tested on Problem B) ---
        function testCostProfileSumMatchesTotalCost(tc)
            % Check that totalCost = sum of stage costs
            prob = tc.solveProbB();
            tc.verifyEqual(sum(prob.CostProfile), prob.totalCost, 'AbsTol', 1e-12)
        end

        function testCostProfileLength(tc)
            % Check that CostProfile has Nstages stage costs + 1 terminal cost
            prob = tc.solveProbB();
            tc.verifyNumElements(prob.CostProfile, prob.Nstages + 1)
        end

        function testStateProfileLength(tc)
            % Check that StateProfile has initial state + Nstages transitions
            prob = tc.solveProbB();
            tc.verifyNumElements(prob.StateProfile{1}, prob.Nstages + 1)
        end

        function testControlProfileLength(tc)
            % Check that ControlProfile has Nstages controls
            prob = tc.solveProbB();
            tc.verifyNumElements(prob.ControlProfile{1}, prob.Nstages)
        end

    end

    methods (Access = private)
        function prob = solveProbA(~)
            % 1 stage: x_next = u, cost = u^2, target [0.9, 1.0]
            prob = DynaProg({0:0.1:1}, {0}, {[0.9, 1.0]}, {0:0.1:1}, 1, ...
                @trivial_sys, 'VFPenalty', 'rift', 'Display', 'off');
            prob = run(prob);
        end

        function prob = solveProbB(~)
            % 2 stages: x_next = x+u, cost = u^2, target {1.0}
            prob = DynaProg({0:0.05:2}, {0}, {[1.0, 1.0]}, {0:0.05:1}, 2, ...
                @trivial_additive_sys, 'VFPenalty', 'rift', 'Display', 'off');
            prob = run(prob);
        end
    end
end
