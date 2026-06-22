classdef tInputValidation < matlab.unittest.TestCase
%tInputValidation  Tests for constructor and property setter validation.
%
%   Verifies that:
%     - Non-cell inputs are accepted and auto-converted
%     - Malformed model functions raise the correct error IDs
%     - Invalid property values raise errors
%     - Coarse state grids trigger warnings

    properties (Constant)
        % Common SV and CV grids, initial state, final state, # stages.
        xGrid = {0:0.1:1} 
        uGrid = {0:0.1:1}
        x0 = {0}
        xf = {[0.9, 1.0]}
        N  = 1
    end

    methods (Test)
        % --- Input auto-conversion ---
        function testNonCellStateGrid_IsAccepted(tc)
            % Test that a scalar state grid is accepted even if it is not
            % cell
            prob = DynaProg(0:0.1:1, tc.x0, tc.xf, tc.uGrid, tc.N, @trivial_sys);
            tc.verifyClass(prob.StateGrid, 'cell')
        end

        function testNonCellStateInitial_IsAccepted(tc)
            % Test that a scalar initial state is accepted even if it is
            % not cell
            prob = DynaProg(tc.xGrid, 0, tc.xf, tc.uGrid, tc.N, @trivial_sys);
            tc.verifyClass(prob.StateInitial, 'cell')
        end

        % --- Model function signature validation ---
        function testModelFun_TwoInputsNoExo_IsAccepted(tc)
            % Test case with two inputs and no exogenous inputs
            fn = @(x, u) trivial_sys(x, u, []);
            prob = DynaProg(tc.xGrid, tc.x0, tc.xf, tc.uGrid, tc.N, fn);
            solved = run(prob);
            tc.verifyNotEmpty(solved.ControlProfile)
        end

        function testModelFun_TwoInputsWithExo_Errors(tc)
            % The case with two inputs and exogenous inputs must error
            fn = @(x, u) trivial_sys(x, u, []);
            prob = DynaProg(tc.xGrid, tc.x0, tc.xf, tc.uGrid, tc.N, fn, ...
                'ExogenousInput', zeros(tc.N, 1));
            tc.verifyError(@() run(prob), 'DynaProg:invalidModelInput')
        end

        function testModelFun_TooFewOutputs(tc)
            % A function with two outputs must error
            prob = DynaProg(tc.xGrid, tc.x0, tc.xf, tc.uGrid, tc.N, @wrong_outputs_sys);
            tc.verifyError(@() run(prob), 'DynaProg:invalidModelOutput')
        end

        % --- StateFinal validation ---
        function testStateFinalWrongOrder(tc)
            % If the final state bounds are not sorted in ascending order,
            % must raise an error
            tc.verifyError( ...
                @() DynaProg(tc.xGrid, tc.x0, {[1.0, 0.5]}, tc.uGrid, tc.N, @trivial_sys), ...
                'DynaProg:wrongOrderStateFinal')
        end

        % --- Grid coarseness warnings ---
        function testTooCoarseStateGrid_WarnsOnRun(tc)
            % Target [0.45, 0.46] has no grid points in grid 0:0.1:1; must
            % raise a warning
            prob = DynaProg(tc.xGrid, tc.x0, {[0.45, 0.46]}, tc.uGrid, tc.N, ...
                @trivial_sys, 'Display', 'warn');
            tc.verifyWarning(@() run(prob), 'DynaProg:tooCoarseStateGrid')
        end
    end
end
