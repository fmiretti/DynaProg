function obj = forward(obj)
% Run the optimization algorithm forward phase

if obj.SafeMode
    control = obj.ControlCombGrid;
    % Vector dimensions corresponding to cvs
    vecdim_cv = 1:length(obj.N_CV);
else
    control = obj.ControlGrid;
    % Vector dimensions corresponding to cvs
    vecdim_cv = (1:length(obj.N_CV)) + length(obj.N_SV);
end

% Progress Bar
if obj.DisplayProgressbar
    fprintf('DP forward progress:    ')
end

% Preallocate profiles
StateProfileMat = zeros(length(obj.StateGrid), obj.Nstages+1);
StateProfileMat(:,1) = [obj.StateInitial{:}];
ControlProfileMat = zeros(length(obj.ControlGrid), obj.Nstages);
obj.CostProfile = zeros(1, obj.Nstages);
% Initialize the state
state = obj.StateInitial;
state_exp = cell(1, length(obj.N_SV));
exoInput = cell(1, size(obj.ExogenousInput, 2));
% Initialize warnings
unfeasFwdWarn = false;
unfeasFwdWarnStages = [];

for k = 1:obj.Nstages
    % Progress Bar
    if obj.DisplayProgressbar
        fprintf('%s%2d %%', ones(1,4)*8, floor((k-1)/obj.Nstages*100));
    end

    % Expand current state to the combined cv grid
    if obj.SafeMode
        for n = 1:length(state)
            state_exp{n} = state{n} + zeros(size(obj.ControlCombGrid{1}));
        end
    else
        state_exp = state;
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
        exoInput_scalar = num2cell(obj.ExogenousInput(k,:));
        if obj.SafeMode
            for n = 1:length(exoInput_scalar)
                exoInput{n} = exoInput_scalar{n}.*ones(size(obj.ControlCombGrid{1}));
            end
        else
            exoInput = exoInput_scalar;
        end
    else
        exoInput = [];
        exoInput_scalar = [];
    end

    % Evaluate state update and stage cost
    [state_next, stageCost, unfeas] = model_wrapper(obj, state_exp, control, exoInput, intVars);
    unfeas = logical(unfeas);
    if obj.UseSplitModel
        unfeas = unfeas | unfeasExt;
    end

    % Expand updated states and unfeas to the combined cv grid
    if ~obj.SafeMode
        for n = 1:length(state_next)
            state_next{n} = state_next{n} + zeros([ones(1, length(obj.N_SV)) obj.N_CV]);
        end
        stageCost = stageCost + zeros([ones(1, length(obj.N_SV)) obj.N_CV]);
        unfeas = unfeas | false([ones(1, length(obj.N_SV)) obj.N_CV]);
    end
    % Enforce state grids
    if obj.EnforceStateGrid
        for n = 1:length(obj.N_SV)
            unfeas(state_next{n} > obj.StateGrid{n}(end) | state_next{n} < obj.StateGrid{n}(1)) = true;
        end
    end
    
    % Find the optimal cvs
    [obj, cv_opt, intVars_opt] = optimalControl(obj, k, state_next, stageCost, unfeas, vecdim_cv, intVars);

    % Advance the simulation
    [state, stageCost_opt, unfeas_opt, addout] = model_wrapper(obj, state, cv_opt, exoInput_scalar, intVars_opt);

    % Check solution validity
    if unfeas_opt
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
    obj.CostProfile(k) = stageCost_opt;
end

% Add terminal cost
obj.CostProfile(end+1) = obj.TerminalCost(state);

% Set cost to infty if the backward phase failed
if obj.failedBackward > 0
    obj.CostProfile = inf(size(obj.CostProfile));
end

% Evaluate total cost
obj.totalCost = sum(obj.CostProfile);

% Forward phase warnings
if obj.DisplayWarnings
    % Check terminal state constraints
    for n = 1:length(obj.StateFinal)
        if ~isempty(obj.StateFinal{n})
            if state{n} < obj.StateFinal{n}(1) || state{n} > obj.StateFinal{n}(2)
                fprintf('\n')
                warning('DynaProg:failedTerminalState', ['The solution violates your terminal state constraints. This may be caused by: \n' ...
                    ' 1) your problem is overconstrained. Try widening the final state constraint bounds.\n' ...
                    ' 2) The state variables grid might be too coarse. Try refining the grids.\n' ...
                    ' 3) You set VFPenalty to ''linear'' but the VFPenFactors are not large enough. Try increasing them.\n' ...
                    'You can also try using the UseLevelSet option.\n' ...
                    ' 4) You set VFPenalty to ''none''. In this case you probably know what you are doing.\n'])
                obj.DisplayProgressbar = false;
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
end

% Store state and control profiles
obj.StateProfile = num2cell(StateProfileMat,2);
obj.ControlProfile = num2cell(ControlProfileMat,2);

% Progress Bar
if obj.DisplayProgressbar
    fprintf('%s%2d %%\n', ones(1,4)*8, 100);
end

end
