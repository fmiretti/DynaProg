function [obj, cv_opt_idx, cost_opt] = updateVF(obj, k, states_next, stageCost, unfeas, vecdim_cv)
%updateVF value function update

% Update Level-Set function
if obj.UseLevelSet
    % Read L(k+1)
    LevelSet_next = obj.LevelSet{k+1}(states_next{:});
    % set LevelSet_next to inf for the unfeasible cvs
    LevelSet_next(unfeas) = obj.myInf;
    % Update level-set function and find L-minimizing cvs
    [LevelSetValue, MinLevelSetCV] = obj.minfun(LevelSet_next, vecdim_cv);
    % Check if the set of reachable cvs U^R(x_k) (cvs that
    % lead to a reachable state) is empty. Note: U^R(x_k) is
    % a subset of U(x_k) (set of feasible cvs).
    isempty_UR = LevelSetValue>0;
    % Level set failure
    if obj.DisplayWarnings && all( LevelSetValue(:) > 0 )
        obj.failedBackward = k;
        obj.DisplayWarnings = false;
        fprintf('\n')
        warning('DynaProg:failedLevelSetUpdate', 'The level set function update has failed: no feasible solution was found. Your problem might be overconstrained or the state variables grid might be too coarse.\n')
    end
    % Construct L approximating function for the current timestep
    obj.LevelSet{k} = griddedInterpolant(obj.StateGridCol, ...
        LevelSetValue, 'linear');
end

% Read VF(k+1)
VF_next =  obj.VF{k+1}(states_next{:});
cost = stageCost + VF_next;
% Set cost-to-go to inf for the unfeasible/unreachable cvs
cost(unfeas) = obj.myInf;
% Find optimal control as a function of the current state
[cost_opt, cv_opt_idx] = obj.minfun(cost, vecdim_cv);

if obj.UseLevelSet
    % For those state grid points where no feasible cv was found
    % (isempty_UR), calculate the VF based on the cv that minimizes
    % the level-set function.
    cost_MinLevelSetCV = cost(MinLevelSetCV);
    cost_opt(isempty_UR) = cost_MinLevelSetCV(isempty_UR);
    if obj.StoreControlMap
        cv_opt_idx(isempty_UR) = MinLevelSetCV(isempty_UR);
    end
end

% Construct VF approximating function for the current timestep
obj.VF{k} = griddedInterpolant(obj.StateGridCol, cost_opt, ...
    'linear');

% Warn the user if there are no feasible trajectories for
% the tail subproblem; do not trigger more than once
if obj.DisplayWarnings && all(cost_opt(:) >= obj.myInf)
    obj.failedBackward = k;
    obj.DisplayWarnings = false;
    fprintf('\n')
    warning('DynaProg:failedCostToGoUpdate', ['The Cost-to-Go update has' ...
        ' failed at stage %d. Your problem might be overconstrained.'], k)
    if obj.DisplayProgressbar
        fprintf('....')
    end
end
end
