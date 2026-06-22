function testModelFun(obj)
%testModelFun Test that the model function returns the correct number of states.
%   The function uses scalar inputs. It does not check the validity of
%   array expansions.
%   Testing split models is also currently unsupported.

if obj.UseSplitModel
    return
end

state = obj.StateInitial;

% First control grid point (scalar)
control = cellfun(@(x) x(1), obj.ControlGrid, 'UniformOutput', false);

% Exogenous input for stage 1 (scalar)
if obj.UseExoInput
    exoInput = num2cell(obj.ExogenousInput(1,:));
else
    exoInput = [];
end

% Call the model and inspect state output
state_next = model_wrapper(obj, state, control, exoInput, []);

if length(state_next) ~= length(obj.N_SV)
    error('DynaProg:wrongSizeStateOutput', ...
        ['You defined %d state variables through StateGrid. ' ...
         'However, the model function returns %d state variables.'], ...
        length(obj.N_SV), length(state_next))
end
end
