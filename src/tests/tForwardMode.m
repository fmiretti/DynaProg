classdef tForwardMode < matlab.unittest.TestCase
%tForwardMode  Tests for ForwardMode and ControlType behaviour.
%
%   Verifies that:
%     - The policyBased algorithm causes ControlMap to be populated after run
%     - valueBased and policyBased give the same totalCost for discrete controls
%     - policyBased with continuous controls gives totalCost within 5% of
%     valueBased

    properties
        ProbValueBasedDiscrete
        ProbPolicyBasedDiscrete
        ProbValueBasedContinuous
        ProbPolicyBasedContinuous
    end

    methods (TestClassSetup)
        function setupProblems(tc)
            grid = {{0:0.1:1}, {0}, {[0.9, 1.0]}, {0:0.1:1}, 1, @trivial_sys, ...
                'VFPenalty', 'rift', 'Display', 'off'};

            tc.ProbValueBasedDiscrete = run(DynaProg(grid{:}, 'ControlType', {'discrete'}));
            tc.ProbPolicyBasedDiscrete = run(DynaProg(grid{:}, 'ControlType', {'discrete'}, ...
                'ForwardMode', 'policyBased'));
            tc.ProbValueBasedContinuous = run(DynaProg(grid{:}, 'ControlType', {'continuous'}));
            tc.ProbPolicyBasedContinuous = run(DynaProg(grid{:}, 'ControlType', {'continuous'}, ...
                'ForwardMode', 'policyBased'));
        end
    end

    methods (Test)
        function testPolicyBased_ControlMapPopulated(tc)
            % After a policyBased run, ControlMap must be non-empty
            tc.verifyFalse(isempty(tc.ProbPolicyBasedDiscrete.ControlMap))
        end

        function testValueBased_ControlMapEmpty(tc)
            % With valueBased (default), ControlMap should stay empty
            tc.verifyTrue(isempty(tc.ProbValueBasedDiscrete.ControlMap))
        end

        function testPolicyBased_SameCostAsValueBased_DiscreteControl(tc)
            % For discrete controls, both modes must produce the same total cost
            tc.verifyEqual(tc.ProbPolicyBasedDiscrete.totalCost, ...
                tc.ProbValueBasedDiscrete.totalCost, 'AbsTol', 1e-6)
        end

        function testPolicyBased_ContinuousControl_SameCostAsValueBased(tc)
            % With continuous interpolation, policyBased should achieve
            % similar cost as ValueBased
            tc.verifyEqual(tc.ProbValueBasedDiscrete.totalCost, ...
                tc.ProbPolicyBasedContinuous.totalCost, 'RelTol', 0.05)
        end

        function testPolicyBased_TerminalStateInBounds(tc)
            % The policyBased forward pass must still satisfy the final state constraint
            xf = tc.ProbPolicyBasedDiscrete.StateProfile{1}(end);
            tc.verifyGreaterThanOrEqual(xf, 0.9)
            tc.verifyLessThanOrEqual(xf, 1.0)
        end
    end
end
