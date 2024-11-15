function [obj, cv_opt, intVars_opt] = optimalControlPolicyBased(obj, k, state, intVars)
%optimalControlPolicyBased find optimal controls for the current stage
%   using control maps

% Get the CVs from the control map
cv_opt =  cellfun(@(x) x( state{:} ), obj.ControlMap(:,k), 'UniformOutput', false);

% Extract the intermediate variables for the optimal cv
if obj.UseSplitModel
    intVars = cellfun(@(x) x .* ones(size(cost)), intVars, 'UniformOutput', false);
    intVars_opt = cellfun(@(x) x(index_opt), intVars, 'UniformOutput', false);
else
    intVars_opt = [];
end

end