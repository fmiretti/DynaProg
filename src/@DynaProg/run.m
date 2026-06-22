function obj = run(obj)
% Check the MATLAB version
if verLessThan('matlab','9.6') % aka 2019a
    obj.minfun =  @(X, vecdim) obj.min_compatibility(X, vecdim);
end
% Validate the model function signature(s). This runs here rather than from
% the property setters so that UseExoInput is final (ExogenousInput is set
% after the system function in the constructor). checkModelFun also sets
% NumAddOutputs, which backward and forward rely on.
if obj.UseSplitModel
    obj = checkModelFun(obj, obj.SysNameExt, 'ext');
    obj = checkModelFun(obj, obj.SysNameInt, 'int');
else
    obj = checkModelFun(obj, obj.SysName, 'single');
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
