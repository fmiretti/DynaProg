function [states_next, stageCost, unfeas, addOutputs] = model_wrapper(obj, state, control, exoInput, IntermediateVars)
getAddOuts = nargout > 3;
if getAddOuts
    addOutputs = cell(1, obj.NumAddOutputs);
end
% Run the function model
if obj.UseSplitModel
    if ~getAddOuts && ~obj.UseExoInput
        [states_next, stageCost, unfeas] = obj.SysNameInt(state, control, IntermediateVars);
        addOutputs = {};
    elseif  ~getAddOuts && obj.UseExoInput
        [states_next, stageCost, unfeas] = obj.SysNameInt(state, control, exoInput, IntermediateVars);
        addOutputs = {};
    elseif getAddOuts && ~obj.UseExoInput
        [states_next, stageCost, unfeas, addOutputs{:}] = obj.SysNameInt(state, control, IntermediateVars);
    else
        [states_next, stageCost, unfeas, addOutputs{:}] = obj.SysNameInt(state, control, exoInput, IntermediateVars);
    end
else
    if ~getAddOuts && ~obj.UseExoInput
        [states_next, stageCost, unfeas] = obj.SysName(state, control);
        addOutputs = {};
    elseif  ~getAddOuts && obj.UseExoInput
        [states_next, stageCost, unfeas] = obj.SysName(state, control, exoInput);
        addOutputs = {};
    elseif  getAddOuts && ~obj.UseExoInput
        [states_next, stageCost, unfeas, addOutputs{:}] = obj.SysName(state, control);
    else
        [states_next, stageCost, unfeas, addOutputs{:}] = obj.SysName(state, control, exoInput);
    end
end
if isempty(unfeas)
    unfeas = zeros(size(state{1}));
end
end
