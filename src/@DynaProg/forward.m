function obj = forward(obj)
% Run the optimization algorithm forward phase

if obj.SafeMode
    control = obj.ControlCombGrid;
    % Vector dimensions corresponding to CVs
    vecdim_cv = 1:length(obj.N_CV);
else
    control = obj.ControlGrid;
    % Vector dimensions corresponding to CVs
    vecdim_cv = (1:length(obj.N_CV)) + length(obj.N_SV);
end

% Progress Bar
progressbar = true;
fprintf('DP forward progress:    ')
% Preallocate profiles
StateProfileMat = zeros(length(obj.StateGrid), obj.Nstages+1);
StateProfileMat(:,1) = [obj.StateInitial{:}];
ControlProfileMat = zeros(length(obj.ControlGrid), obj.Nstages);
% Initialize the state
state = obj.StateInitial;
obj.CostProfile(1) = 0;
state_next = cell(1, length(obj.N_SV));
exoInput = cell(1, size(obj.ExogenousInput, 2));
% Initialize warnings
unfeasFwdWarn = false;
unfeasFwdWarnStages = [];

for k = 1:obj.Nstages
    % Progress Bar
    if progressbar
        fprintf('%s%2d %%', ones(1,4)*8, floor((k-1)/obj.Nstages*100));
    end

    % Expand current state to the full cv grid
    if obj.SafeMode
        for n = 1:length(state)
            state_next{n} = state{n} + zeros(size(obj.ControlCombGrid{1}));
        end
    else
        state_next = state;
    end
    % Create intermediate vaiables
    if obj.UseSplitModel
        intVars = obj.IntermediateVars{k};
        unfeasExt = obj.unFeasExt{k};
    else
        intVars = [];
    end
    % Create exogenous inputs
    if obj.UseExoInput
        exoInput_temp = obj.ExogenousInput(k,:);
        for n = 1:length(exoInput_temp)
            if obj.SafeMode
                exoInput{n} = exoInput_temp(n).*ones(size(obj.ControlCombGrid{1}));
            else
                exoInput{n} = exoInput_temp(n);
            end
        end
    else
        exoInput = [];
    end

    % Evaluate state update and stage cost
    [state_next, stageCost, unFeas] = model_wrapper(obj, state_next, control, exoInput, intVars);
    unFeas = logical(unFeas);
    if obj.UseSplitModel
        unFeas = unFeas | unfeasExt;
    end

    % Expand updated states and unfeas to the full cv grid
    if ~obj.SafeMode
        for n = 1:length(state_next)
            state_next{n} = state_next{n} + zeros([ones(1, length(obj.N_SV)) obj.N_CV]);
        end
        stageCost = stageCost + zeros([ones(1, length(obj.N_SV)) obj.N_CV]);
        unFeas = unFeas | false([ones(1, length(obj.N_SV)) obj.N_CV]);
    end

    % Get Level Set-minimizing CV
    if obj.UseLevelSet
        % Read L(k+1)
        LevelSet_next = obj.LevelSet{k+1}(state_next{:});
        % set LevelSet_next to inf for the unfeasible CVs
        LevelSet_next(unFeas) = obj.myInf;
        % Determine if U^R(x_k) is empty
        isempty_UR = all(LevelSet_next(:) > 0);
        % Find L-minimizing u
        [~, MinLevelSetCV] = obj.minfun(LevelSet_next, vecdim_cv);
    end

    % Read VF(k+1)
    VF_next =  obj.VF{k+1}(state_next{:});
    cost = stageCost + VF_next;
    if obj.UseLevelSet
        cost(unFeas) = obj.myInf;
        cost(LevelSet_next > 0) = obj.myInf;
    end
    % Set cost-to-go to inf for the unfeasible/unreachable CVs
    cost(unFeas) = obj.myInf;

    % Find optimal control as a function of the current state
    [~, index_opt] = obj.minfun(cost, vecdim_cv);
    if obj.UseLevelSet
        % If no reachable cv was found (isempty_UR), select the cv that
        % minimizes the level-set function.
        index_opt(isempty_UR) = MinLevelSetCV;
    end
    cv_opt =  cellfun(@(x) x(index_opt), obj.ControlCombGrid, 'UniformOutput', false);

    % Extract the exogenous inputs for the optimal cv
    if obj.UseExoInput && obj.SafeMode
        exoInput = cellfun(@(x) x(index_opt), exoInput, 'UniformOutput', false);
    end
    % Extract the intermediate variables for the optimal cv
    if obj.UseSplitModel
        intVars = cellfun(@(x) x .* ones(size(cost)), intVars, 'UniformOutput', false);
        intVars_opt = cellfun(@(x) x(index_opt), intVars, 'UniformOutput', false);
    else
        intVars_opt = [];
    end
    % Advance the simulation
    [state, stageCost, unfeas, addout] = model_wrapper(obj, state, cv_opt, exoInput, intVars_opt);

    % Check solution validity
    if unfeas
        unfeasFwdWarnStages(end+1) = k;
        unfeasFwdWarn = true;
    end

    % Update the profiles
    StateProfileMat(:,k+1) = [state{:}]';
    ControlProfileMat(:,k) = [cv_opt{:}]';
    if ~isempty(addout)
        for n = 1:length(addout)
            obj.AddOutputsProfile{n}(k) = addout{n};
        end
    end
    obj.CostProfile(k+1) = stageCost;
end
% Check terminal state constraints
for n = 1:length(obj.StateFinal)
    if ~isempty(obj.StateFinal{n})
        if state{n} < obj.StateFinal{n}(1) || state{n} > obj.StateFinal{n}(2)
            fprintf('\n')
            warning('DynaProg:failedTerminalState', ['The solution violates your terminal state constraints. This may be caused by: \n' ...
                ' 1) your problem is overconstrained. Try widening the final state constraint bounds.\n' ...
                ' 2) The state variables grid might be too coarse. Try refining the grids.\n' ...
                ' 3) You set VFPenalty to ''linear'' but the VFPenFactors are not large enough. Try increasing them.\n' ...
                'You can also try using the Level Set option.\n'])
            progressbar = false;
            break
        end
    end
end
% Print information about constraints violation in the fwd run
if unfeasFwdWarn
    unfeasFwdWarnStages(end+1) = k;
    unfeasFwdWarnStages = unfeasFwdWarnStages(1:min(10, numel(unfeasFwdWarnStages)));
    str = "The solution violates your constraints at stages:\n";
    for n = 1:numel(unfeasFwdWarnStages)
        str = str + "%d ";
    end
    if numel(unfeasFwdWarnStages) >= 10
        str = str + "and others more";
    end
    str = str + ".\nYour problem might be overconstrained or the state variables grid might be too coarse.";
    unfeasFwdWarnStages = num2cell(unfeasFwdWarnStages);
    warning('DynaProg:failedForward', str, unfeasFwdWarnStages{:})
end
% Store state and control profiles
obj.StateProfile = num2cell(StateProfileMat,2);
obj.ControlProfile = num2cell(ControlProfileMat,2);
% Progress Bar
if progressbar
    fprintf('%s%2d %%\n', ones(1,4)*8, 100);
end
end
