function obj = backward(obj)
% Run the optimization algorithm backward phase
if obj.SafeMode
    state = obj.StateFullGrid;
    control = obj.ControlFullGrid;
else
    state = obj.StateGrid;
    control = obj.ControlGrid;
end

% Progress Bar 
if obj.DisplayProgressbar
    fprintf('DP backward progress:    ')
end
% Vector dimensions corresponding to CVs
vecdim_cv = (length(obj.N_SV)+1):(length(obj.N_CV)+length(obj.N_SV));

% Backward Loop
for k = obj.Nstages:-1:1
    % Progress Bar
    if obj.DisplayProgressbar
        fprintf('%s%2d %%', ones(1,4)*8, floor((obj.Nstages-k)/obj.Nstages*100));
    end

    if obj.UseSplitModel
        if obj.SafeMode
            % Expand the intermediate variables to the combined grid
            intVars = cellfun(@(x) repmat(shiftdim(x, -length(obj.N_SV)), [obj.N_SV ones(size(obj.N_CV))]), obj.IntermediateVars{k}, 'UniformOutput', false);
            unfeasExt = repmat(shiftdim(obj.unFeasExt{k}, -length(obj.N_SV)), [obj.N_SV ones(size(obj.N_CV))]);
        else
            intVars = obj.IntermediateVars{k};
            unfeasExt = obj.unFeasExt{k};
        end
    else
        intVars = [];
    end
    % Create exogenous inputs
    if obj.UseExoInput
        currentExoInput = obj.ExogenousInput(k,:);
        for n = 1:length(currentExoInput)
            if obj.SafeMode
                exoInput{n} = currentExoInput(n).*ones(size(obj.StateFullGrid{1}));
            else
                exoInput{n} = currentExoInput(n);
            end
        end
    else
        exoInput = [];
    end

    % State update
    [states_next, stageCost, unFeasInt] = model_wrapper(obj, state, control, exoInput, intVars);
    unFeasInt = logical(unFeasInt);

    % Model output checks
    if k == obj.Nstages
        if iscell(stageCost)
            error('DynaProg:wrongFormatStageCost', ['The stage '...
                'cost must be returned as a numeric type, not a cell.'])
        elseif ~isnumeric(stageCost)
            error('DynaProg:wrongFormatStageCost', ['The stage '...
                'cost must be returned as a numeric type.'])
        end
    end

    if obj.DisplayWarnings && all(unFeasInt(:) == true)
        obj.failedBackward = k;
        obj.DisplayWarnings = false;
        fprintf('\n');
        warning('DynaProg:unfeasModel', ['There are no feasible state/controls '...
            'at stage %d. Check the model function, state grids and control grids.'], k)
        if obj.DisplayProgressbar
            fprintf('....')
        end
    end

    % feasibility: include external model unfeasibility
    if obj.UseSplitModel
        unfeas = unFeasInt | unfeasExt;
    else
        unfeas = unFeasInt;
    end
    if ~obj.SafeMode
        % Expand updated states and unfeas to the combined grid
        for n = 1:length(states_next)
            states_next{n} = states_next{n} + zeros([obj.N_SV obj.N_CV]);
        end
        stageCost = stageCost + zeros([obj.N_SV obj.N_CV]);
        unfeas = unfeas | false([obj.N_SV obj.N_CV]);
    end
    % Enforce state grids
    if obj.EnforceStateGrid
        for n = 1:length(obj.N_SV)
            unfeas(states_next{n} > obj.StateGrid{n}(end) | states_next{n} < obj.StateGrid{n}(1)) = true;
        end
    end

    % Update the value function
    [obj, cv_opt] = updateVF(obj, k, states_next, stageCost, unfeas, vecdim_cv);

    % Store cv map
    if obj.StoreControlMap
        for n = 1:length(obj.N_CV)
            cv_opt_sub = cell(1, length(obj.N_SV) + length(obj.N_CV));
            [cv_opt_sub{:}] = ind2sub([obj.N_SV, obj.N_CV], cv_opt);
            obj.ControlMap{n,k} = squeeze(obj.ControlGrid{n}(cv_opt_sub{n+length(obj.N_SV)}));
        end
    end

end

if obj.DisplayWarnings && ( obj.VF{1}(obj.StateInitial) > obj.myInf )
    obj.failedBackward = 1;
    obj.DisplayWarnings = false;
    fprintf('\n')
    warning('DynaProg:unfeasibleInitialState', ['There is no feasible trajectory ' ...
        'starting from the initial state. Check the initial state. ' ...
        'Also, your problem might be overcostrained.\n'])
end

% Progress Bar
if obj.DisplayProgressbar
    fprintf('%s%2d %%\n', ones(1,4)*8, 100);
end

end
