function obj = run(obj)
% Check the MATLAB version
if verLessThan('matlab','9.6') % aka 2019a
    obj.minfun =  @(X, vecdim) obj.min_compatibility(X, vecdim);
end
% Create computational grids
obj = create_grids(obj);
if obj.UseSplitModel
    obj = create_intVars(obj);
end
% Generate value functions
obj = backward(obj);
% Generate optimal trajectories
obj = forward(obj);
% Remove VF if its storage was not required
if ~obj.StoreValueFunction
    obj.VF = [];
end
% Reset verbosity flags
obj.Display = obj.Display;

end
