function [states_next, stageCost, unFeas, addOutputs] = model_wrapper(obj, state, control, exoInput, IntermediateVars)
getAddOuts = nargout > 3;
if getAddOuts
    addOutputs = cell(1, obj.NumAddOutputs);
end
% Run the function model
if obj.UseSplitModel
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
if isempty(unFeas)
    unFeas = zeros(size(state{1}));
end
end
