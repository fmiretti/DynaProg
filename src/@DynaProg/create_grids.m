function obj = create_grids(obj)
% Create computational grids, initialize terminal value
% function and level-set function

% checks
if ~isempty(obj.StateGrid)
    obj = check_StateFinal(obj);
    if length(obj.StateFinal) ~= length(obj.StateGrid)
        error('DynaProg:wrongSizeStateFinal', ['You must specify '...
            'one final condition for each of the state variables.']);
    end
end
if length(obj.StateInitial) ~= length(obj.StateGrid)
    error('DynaProg:wrongSizeStateInit', ['You must specify one '...
        'initial condition for each of the state variables.']);
end
% Set the VF initialization method if unspecified
if strcmp(obj.VFInitialization, 'auto')
    if isempty(obj.TerminalCost)
        if obj.UseLevelSet
            obj.VFInitialization = 'linear';
        else
            obj.VFInitialization = 'rift';
        end
    else
        obj.VFInitialization = 'manual';
    end
end
% Checks on the VF initialization method
if ~isempty(obj.TerminalCost) && ~strcmp(obj.VFInitialization, 'manual')
    warning('DynaProg:ignoredTerminalCost', 'You specified a terminal cost with TerminalCost but VFInitialization is not set to ''manual''. Ignoring your terminal cost.' )
end
if isempty(obj.TerminalCost) && strcmp(obj.VFInitialization, 'manual')
    warning('DynaProg:emptyTerminalCost', 'You set VFInitialization to ''manual'' but you did not specify a terminal cost. Set VFInitialization to a valid string or specify a terminal cost with TerminalCost.' )
end
% Set the Level Set initialization method if unspecified
if isempty(obj.LevelSetInitialization)
    obj.LevelSetInitialization = 'linear';
end

% State grids as N_SV_n-by-1 vectors, needed to create VF interpolants
obj.StateGridCol = cellfun(@(x) x(:), obj.StateGrid, 'UniformOutput', false);

% Full SV and CV grids
if obj.SafeMode
    obj.StateFullGrid = cell(1, length(obj.N_SV));
    obj.ControlFullGrid = cell(1, length(obj.N_CV));
    [obj.StateFullGrid{:}, obj.ControlFullGrid{:}] = ndgrid(obj.StateGrid{:}, obj.ControlGrid{:});
end

% Combined CV grid
obj.ControlCombGrid = cell(1, length(obj.ControlGrid));
[obj.ControlCombGrid{:}] = ndgrid(obj.ControlGrid{:});
for n = 1:length(obj.N_CV)
    obj.ControlGrid{n} = obj.ControlGrid{n}(:);
    obj.ControlGrid{n} = shiftdim(obj.ControlGrid{n}, -length(obj.StateGrid) - (n-1));
end
% Initialize terminal VF
StateFullGrid = cell(1, length(obj.StateGrid));
[StateFullGrid{:}] = ndgrid(obj.StateGrid{:});
VFN = zeros(size(StateFullGrid{1}));
if ~isempty(obj.StateFinal)
    switch obj.VFInitialization
        case 'linear' %VFN is proportional to the distance from the target set
            for n = 1:length(obj.StateGrid)
                if ~isempty(obj.StateFinal{n})
                    VFN = VFN + obj.VFFactors(n) .* ( max(obj.StateFinal{n}(1)-StateCombGrid{n}, 0) + max(StateCombGrid{n}-obj.StateFinal{n}(2), 0) );
                end
            end
            VFN(isinf(VFN)) = 0;
        case 'rift' % VFN is inf outside the target set
            for n=1:length(obj.StateGrid)
                if ~isempty(obj.StateFinal{n})
                    VFN( StateCombGrid{n} > obj.StateFinal{n}(2) ...
                        | StateCombGrid{n} < obj.StateFinal{n}(1)) = obj.myInf;
                end
            end
        case 'manual'
            VFN = obj.TerminalCost;
            % Check user-supplied terminal cost
            if ~isequal(size(VFN), obj.N_SV)
                error('DynaProg:wrongSizeTerminalCost', strjoin({'The terminal cost you provided has wrong size. It should be', sprintf('%dx', obj.N_SV), sprintf('\b\b (the lengths of the state grids).'), '\n'}))
            end
    end
end
% Allocate cell array for the VF interpolants and create the interpolant
obj.VF = cell(1, obj.Nstages+1);
obj.VF{end} = griddedInterpolant(obj.StateGridCol, VFN, ...
    'linear');
% Initialize Level Set function
if obj.UseLevelSet
    LevelSetN = cell(1,length(obj.StateGrid));
    for n=1:length(obj.StateGrid)
        if ~isempty(obj.StateFinal{n})
            LevelSetN{n} = max( obj.StateFinal{n}(1) - StateCombGrid{n}, ...
                StateCombGrid{n} - obj.StateFinal{n}(2) );
        else
            LevelSetN{n} = zeros(size(StateCombGrid{n}));
        end
    end
    LevelSetN = cat(length(LevelSetN)+1,LevelSetN{:});
    LevelSetN = max(LevelSetN,[],ndims(LevelSetN));

    % Allocate cell array for the LevelSet interpolants and create the interpolant
    obj.LevelSet = cell(1, obj.Nstages+1);
    obj.LevelSet{end} = griddedInterpolant(obj.StateGridCol, LevelSetN, ...
        'linear');

end
end