function obj = create_intVars(obj)
% Run external model to create intermediate variables
exoInput = cell(1, size(obj.ExogenousInput, 2));
if obj.SafeMode
    control = obj.ControlCombGrid;
else
    control = obj.ControlGrid;
end


for k = 1:obj.Nstages
    % Create exogenous inputs
    if obj.UseExoInput
        currentExoInput = obj.ExogenousInput(k,:);
        for n = 1:length(currentExoInput)
            if obj.SafeMode
                exoInput{n} = currentExoInput(n).*ones(size(obj.ControlCombGrid{1}));
            else
                exoInput{n} = currentExoInput(n);
            end
        end
    else
        exoInput = [];
    end
    % Create intVars for the current timestep
    [obj.IntermediateVars{k}, obj.unFeasExt{k}] = obj.SysNameExt(control, exoInput);
    % Check intermediate variables size
    sizes = cellfun(@(x) size(x), obj.IntermediateVars{k}, 'UniformOutput', false);
    if obj.SafeMode
        wrong_size = ~cellfun(@(x) isequal(x, obj.N_CV), sizes);
    else
        for n = 1:length(sizes)
            sizes_sv = sizes{n}(1:length(obj.N_SV));
            sizes_cv = ones(1, length(obj.N_CV));
            sizes_cv(1:(length(sizes{n})-length(obj.N_SV))) = sizes{n}(length(obj.N_SV)+1:end);
            wrong_size_cv = ~all(sizes_cv == obj.N_CV | sizes_cv == ones(size(obj.N_CV)));
            wrong_size(n) = wrong_size_cv | ~isequal(sizes_sv, ones(size(obj.N_SV)));
        end
    end
    if any(wrong_size)
        error('DynaProg:wrongSizeIntVars', ['Intermediate variables #' num2str(wrong_size) ' in the external model function has a wrong size.'])
    end
end
end