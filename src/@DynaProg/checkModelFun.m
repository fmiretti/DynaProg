function obj = checkModelFun(obj, name, mode)
% Extract the model function name
info = functions(name);
if strcmp(info.type, 'anonymous')
    func_name = regexp(info.function, '\)\w*\(', 'match');
    func_name = func_name{1}(2:end-1);
else
    func_name = name;
end
% Check the model function name
try
    nargout(func_name);
catch ME
    switch ME.identifier
        case {'MATLAB:narginout:functionDoesnotExist', 'MATLAB:narginout:notValidMfile'}
            error(['Model function "' char(func_name) '" not found. Check the function name. Check the filename and path if the function is in a file.'])
        otherwise
            rethrow(ME)
    end
end

% Check the number of inputs
switch mode
    case 'single'
        if nargin(name) < 3
            error('DynaProg:invalidModelInput', ['The model function must have at least three inputs (x, u and w).\n' ...
                'If you do not want to use exogenous inputs, suppress the third input (replace it with a tilde ~)'])
        end
    case 'ext'
        if nargin(name) < 2
            error('DynaProg:invalidModelInput', ['The external model function must have at least two inputs (u and w).\n' ...
                'If you do not want to use exogenous inputs, suppress the second input (replace it with a tilde ~)'])
        end
    case 'int'
        if nargin(name) < 4
            error('DynaProg:invalidModelInput', ['The internal model function must have at least four inputs (x, u, w and intVar).\n' ...
                'If you do not want to use exogenous inputs, suppress the third input (replace it with a tilde ~)'])
        end
end

% Determine the number of additional outputs declared in the
% user's system function signature
switch mode
    case {'single', 'int'}
        if nargout(func_name) < 3
            error('DynaProg:invalidModelOutput', 'The function model has a wrong number of outputs.\n')
        end
        obj.NumAddOutputs = nargout(func_name) - 3;
end
end
