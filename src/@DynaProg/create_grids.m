function obj = create_grids(obj)
% Create computational grids, initialize terminal value
% function and level-set function

% checks
if ~isempty(obj.StateGrid)
    
    if obj.DisplayWarnings
        obj = check_StateFinal(obj);
    end

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
if strcmp(obj.VFPenalty, 'auto')
    if obj.UseLevelSet
        obj.VFPenalty = 'linear';
    else
        obj.VFPenalty = 'rift';
    end
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

% Initialize terminal VF with the terminal cost
VFN = obj.TerminalCost(obj.StateGrid);
% Check user-supplied terminal cost
if ( length(obj.N_SV) > 1 && ~isequal(size(VFN), obj.N_SV) ) || ( length(obj.N_SV) == 1 && ~isequal(length(VFN), obj.N_SV) )
    error('DynaProg:wrongSizeTerminalCost', strjoin({'The terminal cost function you provided returns a cost with wrong size. It should be', sprintf('%dx', obj.N_SV), sprintf('\b\b (the lengths of the state grids).'), '\n'}))
end

% Combined SV grid
StateCombGrid = cell(1, length(obj.StateGrid));
[StateCombGrid{:}] = ndgrid(obj.StateGrid{:});

% Add penalty term to the terminal VF
if ~isempty(obj.StateFinal)
    switch obj.VFPenalty
        case 'linear' %VFN is proportional to the distance from the target set
            for n = 1:length(obj.StateGrid)
                if ~isempty(obj.StateFinal{n})
                    VFN = VFN + obj.VFPenFactors(n) .* ( max(obj.StateFinal{n}(1)-StateCombGrid{n}, 0) + max(StateCombGrid{n}-obj.StateFinal{n}(2), 0) );
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