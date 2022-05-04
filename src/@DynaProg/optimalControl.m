function [obj, cv_opt, intVars_opt] = optimalControl(obj, k, state_next, stageCost, unfeas, vecdim_cv, intVars)
%optimaControl find optimal controls for the current stage

% Get level set - minimizing cv
if obj.UseLevelSet
    % Read L(k+1)
    LevelSet_next = obj.LevelSet{k+1}(state_next{:});
    % Set LevelSet_next to inf for the unfeasible cvs
    LevelSet_next(unfeas) = obj.myInf;
    % Determine if U^R(x_k) is empty
    isempty_UR = all(LevelSet_next(:) > 0);
end

% Read VF(k+1)
VF_next =  obj.VF{k+1}(state_next{:});
cost = stageCost + VF_next;
% Set cost-to-go to inf for the unfeasible cvs
cost(unfeas) = obj.myInf;
% Set the cost to go to inf for cvs leading to unreachable states
if obj.UseLevelSet
    cost(LevelSet_next > 0) = obj.myInf;
end

% Find optimal control as a function of the current state
[~, index_opt] = obj.minfun(cost, vecdim_cv);
if obj.UseLevelSet
    % If no reachable cv was found (isempty_UR), then use the L-minimizing u
    [~, MinLevelSetCV] = obj.minfun(LevelSet_next, vecdim_cv);
    index_opt(isempty_UR) = MinLevelSetCV;
end
cv_opt =  cellfun(@(x) x(index_opt), obj.ControlCombGrid, 'UniformOutput', false);

% Extract the intermediate variables for the optimal cv
if obj.UseSplitModel
    intVars = cellfun(@(x) x .* ones(size(cost)), intVars, 'UniformOutput', false);
    intVars_opt = cellfun(@(x) x(index_opt), intVars, 'UniformOutput', false);
else
    intVars_opt = [];
end

end