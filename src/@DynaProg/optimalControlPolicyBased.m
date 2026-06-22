function [obj, cv_opt, intVars_opt] = optimalControlPolicyBased(obj, k, state, ~)
%optimalControlPolicyBased find optimal controls for the current stage
%   using control maps

% Get the CVs from the control map
cv_opt = cellfun(@(x) x(state{:}), obj.ControlMap(:,k), 'UniformOutput', false);

% Extract the intermediate variables at the optimal cv.
% SysNameExt is called at cv_opt directly so that off-grid values produced
% by linear interpolation of a continuous ControlMap are handled correctly.
if obj.UseSplitModel
    if obj.UseExoInput
        exoInput = num2cell(obj.ExogenousInput(k,:));
    else
        exoInput = [];
    end
    [intVars_opt, ~] = obj.SysNameExt(cv_opt, exoInput);
else
    intVars_opt = [];
end

end