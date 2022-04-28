function [states_next, stageCost, unFeas, addOutputs] = model_wrapper(obj, state, control, exoInput, IntermediateVars)

getAddOuts = nargout > 3;
if getAddOuts
    addOutputs = cell(1, obj.NumAddOutputs);
end

%% Run the function model
if obj.UseSplitModel
    % Run the internal function model
    if ~getAddOuts && ~obj.UseExoInput
        [states_next, stageCost, unFeas] = obj.SysNameInt(state, control, IntermediateVars);
        addOutputs = {};
    elseif  ~getAddOuts && obj.UseExoInput
        [states_next, stageCost, unFeas] = obj.SysNameInt(state, control, exoInput, IntermediateVars);
        addOutputs = {};
    elseif getAddOuts && ~obj.UseExoInput
        [states_next, stageCost, unFeas, addOutputs{:}] = obj.SysNameInt(state, control, IntermediateVars);
    else
        [states_next, stageCost, unFeas, addOutputs{:}] = obj.SysNameInt(state, control, exoInput, IntermediateVars);
    end
else
    % Run the function model
    if ~getAddOuts && ~obj.UseExoInput
        [states_next, stageCost, unFeas] = obj.SysName(state, control);
        addOutputs = {};
    elseif  ~getAddOuts && obj.UseExoInput
        [states_next, stageCost, unFeas] = obj.SysName(state, control, exoInput);
        addOutputs = {};
    elseif  getAddOuts && ~obj.UseExoInput
        [states_next, stageCost, unFeas, addOutputs{:}] = obj.SysName(state, control);
    else
        [states_next, stageCost, unFeas, addOutputs{:}] = obj.SysName(state, control, exoInput);
    end
end

% Set unfeas for unconstrained problems
if isempty(unFeas)
    unFeas = zeros(size(state{1}));
end

%% Model output checks
if k == obj.Nstages
    if iscell(stageCost)
        error('DynaProg:wrongFormatStageCost', ['The stage '...
            'cost must be returned as a numeric type, not a cell.'])
    elseif ~isnumeric(stageCost)
        error('DynaProg:wrongFormatStageCost', ['The stage '...
            'cost must be returned as a numeric type.'])
    end
end

end
