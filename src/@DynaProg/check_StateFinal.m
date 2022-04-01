function obj = check_StateFinal(obj)
for n = 1:length(obj.StateGrid)
    if ~isempty(obj.StateFinal{n})
        if obj.UseLevelSet && (obj.StateFinal{n}(1) == obj.StateFinal{n}(2))
            warning('DynaProg:levelSetFail', 'When the level-set method is enabled, you should not specify the target state as a single point. Make sure you specify distinct lower and upper final state bounds.')
        else
            count = sum(obj.StateGrid{n} >= obj.StateFinal{n}(1) & obj.StateGrid{n} <= obj.StateFinal{n}(2));
            if count == 0
                warning('DynaProg:tooCoarseStateGrid', ['The state variable grid for SV #' num2str(n) ' is too coarse. At least one of the grid points should lie inside the final state bounds.'])
            elseif count < 3
                warning('DynaProg:coarseStateGrid', ['The state variable grid for SV #' num2str(n) ' is very coarse. Consider defining the grid so that at least three of the grid points lie inside the final state bounds, or loosen the final state bounds.'])
            end
        end
    end
end
end