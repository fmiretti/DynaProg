classdef tRegression < matlab.unittest.TestCase
%tRegression  End-to-end regression tests on real example problems.
%
%   Runs the two_tanks and cart (double integrator) examples and checks:
%     - Solution is finite
%     - Terminal state satisfies the specified constraints
%     - Profile structure matches expected lengths
%     - Total cost matches a locked-in reference value
%
%   The reference total cost values were recorded from a known-good run
%   and must be updated if the algorithm changes intentionally.

    properties (Constant)
        % Reference values
        TwoTanks_ExpectedCost = 1.390000000000001;
        TwoTanks_ExpectedFinalState = [0.500044805354877, 0.500044805354877];
        Cart_ExpectedCost = 5.672000000000001;
        Cart_ExpectedFinalState = [0.692, 0.020];
    end

    properties
        ProbTwoTanks
        ProbCart
    end

    methods (TestClassSetup)

        function setupProblems(tc)
            x_grid  = {0:0.02:1, 0:0.02:1};
            x_init  = {0, 0};
            x_final = {[0.5, 1], [0.5, 1]};
            u_grid  = {0:0.05:1, 0:0.05:1};
            tc.ProbTwoTanks = run(DynaProg(x_grid, x_init, x_final, u_grid, 200, ...
                @two_tanks, 'UseLevelSet', true, 'Display', 'off'));

            x_grid  = {0:0.01:1, -0.2:0.01:1.2};
            x_init  = {0, 0};
            x_final = {[0.69, 0.71], [-0.02, 0.02]};
            u_grid  = {-5:0.05:5};
            tc.ProbCart = run(DynaProg(x_grid, x_init, x_final, u_grid, 10, ...
                @cart, 'Display', 'off'));
        end

    end

    methods (Test)
        % --- two_tanks (UseLevelSet = true) ---
        function twoTanks_CostIsFinite(tc)
            tc.verifyTrue(isfinite(tc.ProbTwoTanks.totalCost))
        end

        function twoTanks_TerminalState_SV1_InBounds(tc)
            xf1 = tc.ProbTwoTanks.StateProfile{1}(end);
            tc.verifyGreaterThanOrEqual(xf1, 0.5)
            tc.verifyLessThanOrEqual(xf1, 1.0)
        end

        function twoTanks_TerminalState_SV2_InBounds(tc)
            xf2 = tc.ProbTwoTanks.StateProfile{2}(end);
            tc.verifyGreaterThanOrEqual(xf2, 0.5)
            tc.verifyLessThanOrEqual(xf2, 1.0)
        end

        function twoTanks_StateProfileLength(tc)
            tc.verifyNumElements(tc.ProbTwoTanks.StateProfile{1}, 200 + 1)
        end

        function twoTanks_ControlProfileLength(tc)
            tc.verifyNumElements(tc.ProbTwoTanks.ControlProfile{1}, 200)
        end

        function twoTanks_CostProfileSumMatchesTotalCost(tc)
            tc.verifyEqual(sum(tc.ProbTwoTanks.CostProfile), tc.ProbTwoTanks.totalCost, ...
                'AbsTol', 1e-10)
        end

        function twoTanks_TotalCostMatchesReferecne(tc)
            tc.verifyEqual(tc.ProbTwoTanks.totalCost, tc.TwoTanks_ExpectedCost, 'AbsTol', 1e-10)
        end

        function twoTanks_FinalStateMatchesReference(tc)
            xf = [tc.ProbTwoTanks.StateProfile{1}(end), tc.ProbTwoTanks.StateProfile{2}(end)];
            tc.verifyEqual(xf, tc.TwoTanks_ExpectedFinalState, 'AbsTol', 1e-10)
        end

        % --- cart (double integrator, no exo input) ---
        function cart_CostIsFinite(tc)
            tc.verifyTrue(isfinite(tc.ProbCart.totalCost))
        end

        function cart_TerminalPosition_InBounds(tc)
            xf1 = tc.ProbCart.StateProfile{1}(end);
            tc.verifyGreaterThanOrEqual(xf1, 0.69)
            tc.verifyLessThanOrEqual(xf1, 0.71)
        end

        function cart_TerminalVelocity_InBounds(tc)
            xf2 = tc.ProbCart.StateProfile{2}(end);
            tc.verifyGreaterThanOrEqual(xf2, -0.02)
            tc.verifyLessThanOrEqual(xf2, 0.02)
        end

        function cart_StateProfileLength(tc)
            tc.verifyNumElements(tc.ProbCart.StateProfile{1}, 10 + 1)
        end

        function cart_ControlProfileLength(tc)
            tc.verifyNumElements(tc.ProbCart.ControlProfile{1}, 10)
        end

        function cart_CostProfileSumMatchesTotalCost(tc)
            tc.verifyEqual(sum(tc.ProbCart.CostProfile), tc.ProbCart.totalCost, 'AbsTol', 1e-10)
        end

        function cart_TotalCostMatchesReference(tc)
            tc.verifyEqual(tc.ProbCart.totalCost, tc.Cart_ExpectedCost, 'AbsTol', 1e-10)
        end

        function cart_FinalStateMatchesReference(tc)
            xf = [tc.ProbCart.StateProfile{1}(end), tc.ProbCart.StateProfile{2}(end)];
            tc.verifyEqual(xf, tc.Cart_ExpectedFinalState, 'AbsTol', 1e-10)
        end
    end
end
