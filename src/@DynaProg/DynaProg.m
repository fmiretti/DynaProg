classdef (CaseInsensitiveProperties=true) DynaProg
    %DynaProg  Create a DynaProg problem structure.
    %
    %   prob = DynaProg(STATEGRID, STATEINITIAL, STATEFINAL,
    %   CONTROLGRID, NSTAGES, SYSNAME) creates the basic problem structure prob.
    %
    %   SYSNAME is a function handle to the model function.
    %   The model function must return the updated state value
    %   and the stage cost as a function of the current state value and the
    %   control variables. Additionally, it can accept an exogenous input
    %   and it can return unfeasibilities.
    %   The structure of the model function must be:
    %
    %      [X_NEXT, STAGE_COST, UNFEAS] = SYSNAME(X, U, W);
    %
    %      where:
    %       X is a cell array, where each cell contains the value for a
    %       state variable.
    %       U is a cell array, where each cell contains the value for a
    %       control variable.
    %       W (optional) is a cell array, where each cell contains the
    %       value for an exogenous input variable. If not needed, replace
    %       with a tilde (~) in the function signature.
    %       X_NEXT is a cell array, where each cell contains the value for
    %       the updated state variable.
    %       STAGE_COST is a numeric variable, containing the stage cost.
    %       UNFEAS is either an empty array a boolean variable, set to true
    %       for unfeasible (not allowed) combinations of state variables,
    %       control variables and exogenous inputs. If not needed, return a
    %       empty array ([]) in its place.
    %
    %   STATEGRID is a cell array of numeric vectors, where each vector
    %   defines the discretized grid for a state variable.
    %
    %   STATEINITIAL is a cell array of scalar values, where each value
    %   defines the initial value for a state variable.
    %
    %   STATEFINAL is a cell array of two-element vectors, where each
    %   vector defines lower and upper constraints for the final value for
    %   each state variable. Each vector must contain the lower and upper
    %   bound, in this order.
    %
    %   CONTROLGRID is a cell array of numeric vectors, where each vector
    %   defines the discretized grid for a control variable.
    %
    %   NSTAGES is the number of stages of the optimization problem.
    %
    %   prob = DynaProg(STATEGRID, STATEINITIAL, [],
    %   CONTROLGRID, NSTAGES, SYSNAME) creates the basic problem structure
    %   prob, with no contraints on the final state.
    %
    %   prob = DynaProg(STATEGRID, STATEINITIAL, STATEFINAL,
    %   CONTROLGRID, NSTAGES, EXTSYSNAME, INTSYSNAME) creates the basic
    %   problem structure prob, using the split system approach.
    %
    %   EXTSYSNAME is a function handle to the external model
    %   function. The external model function must return the intermediate
    %   variables as a function of the control variables.
    %   Additionally, it can accept exogenous inputs and it can return
    %   unfeasibilities.
    %   The structure of the external model function must be:
    %
    %      [M, UNFEAS] = EXTSYSNAME(U, W);
    %
    %      where:
    %       U is a cell array, where each cell contains the value for a
    %       control variable.
    %       W (optional) is a cell array, where each cell contains the
    %       value for an exogenous input variable.
    %       M is a cell array, where each cell contains the value for
    %       an intermediate variable.
    %       UNFEAS (optional) is a boolean variable, set to true for
    %       unfeasible (not allowed) combinations of state variables,
    %       control variables and exogenous inputs. If not needed, return a
    %       empty array ([]) in its place.
    %
    %   INTSYSNAME is a function handle to the internal model function.
    %   The internal model function must return the updated state value
    %   and the stage cost as a function of the current state value, the
    %   control variables and the intermediate variables.
    %   Additionally, it can accept an exogenous input and it can return
    %   unfeasibilities (both are optional).
    %   The structure of the model function must be:
    %
    %      [X_NEXT, STAGE_COST, UNFEAS] = INTSYSNAME(X, U, M, W);
    %
    %      where:
    %       X is a cell array, where each cell contains the value for a
    %       state variable.
    %       U is a cell array, where each cell contains the value for a
    %       control variable.
    %       M is a cell array, where each cell contains the value for
    %       an intermediate variable.
    %       W (optional) is a cell array, where each cell contains the
    %       value for an exogenous input variable.
    %       X_NEXT is a cell array, where each cell contains the value for
    %       the updated state variable.
    %       STAGE_COST is a numeric variable, containing the stage cost.
    %       UNFEAS (optional) is a boolean variable, set to true for
    %       unfeasible (not allowed) combinations of state variables,
    %       control variables and exogenous inputs. If not needed, return a
    %       empty array ([]) in its place.
    %
    %   prob = DynaProg(_, 'PROPERTY1', VALUE, ..., 'PROPERTYN', VALUE)
    %   specifies additional properties and information with
    %   parameter/value pairs.
    %
    %   Valid properties are:
    %       UseLevelSet: enable Level-Set DP. Defaults to false if
    %           unspecified.
    %       ExogenousInput: specify exogenous inputs required for your
    %           model in a cell array of numeric vectors. Each vector must have
    %           the same length as the number of stages of the optimization
    %           problem.
    %       SafeMode: enables Safe Mode (see the documentation)
    %       myInf: if VFPenalty is set to 'rift', myInf defines the
    %           penalty cost for unfeasible terminal states.
    %       EnforceStateGrid: set a constraint so that the state variables
    %           do not exceed the state grids. Defaults to true.
    %       Display: level of display in the command window.
    %           - 'off' displays no output.
    %           - 'warn' displays only warnings and hides the progress bar.
    %           - 'detailed' (default) displays both warnings and the 
    %             progress bar.
    %     # Terminal Cost (Value Function initialization) settings
    %       TerminalCost: define a custom terminal cost as a function 
    %           handle. The function handle must accept only one input 
    %           (the final states in a cell array) and it must return the 
    %           terminal cost. The size of the output must be composed by 
    %           the lenghts of the state variable grids (i.e.
    %           length(x1_grid) x length(x2_grid) x ... x length(xn_grid)).
    %       VFPenalty: specify how final state values outside of
    %           the final state constraints bounds should be penalized. Set
    %           to 'rift' to penalize them with a myInf term. Set to
    %           'linear' to penalize them with a term proportional to the
    %           distance from the bounds. Set to 'none' to disable the
    %           penalty altogether.
    %       VFPenFactors: if VFPenalty is set to 'linear', VFPenFactors
    %           define the proportionality factor for each sv. Specify as a
    %           numeric array.
    %     # Plots settings
    %       StateName: specify state variables names in a string
    %           array.
    %       ControlName: specify state variables names in a string
    %           array.
    %       CostName: specify the cumulative cost name as a string.
    %       Time: specify time instead of stages. This property is only
    %           used in the plots produced with the plot method. It does not
    %           affect the optimization.
    %
    %   Author: Federico Miretti
    %
    % <a href="matlab:web html/index.html">Open the full documentation</a>

    %% Properties
    properties
        % Constructor: Main properties
        SysName function_handle
        SysNameExt function_handle
        SysNameInt function_handle
        StateGrid
        StateInitial
        StateFinal = [];
        ControlGrid
        Nstages
        % Constructor: Name-Value pair arguments
        ExogenousInput = [];
        UseLevelSet logical = false;
        StoreValueFunction logical = false;
        StoreControlMap logical = false;
        SafeMode = false;
        Time double = [];
        myInf double = 1e6;
        EnforceStateGrid logical = true;
        VFPenalty char {mustBeMember(VFPenalty, {'none', 'rift', 'linear', 'auto'})} = 'auto';
        LevelSetInitialization char = [];
        TerminalCost function_handle = @(x) zeros( [cellfun(@length, x) 1] );
        % Results
        totalCost
        StateProfile
        ControlProfile
        CostProfile
        AddOutputsProfile
        ControlMap
    end
    properties (Access = protected)
        % Computational grids, value and level-set function
        StateFullGrid % Full grids for the SVs 
        StateGridCol % State grids as column vectors, needed to create VF interpolants
        ControlFullGrid % Full grids for the CVs
        ControlCombGrid % Combined CV grid (only combines the CV dimensions)
        LevelSet % Level-Set function
        IntermediateVars
        unFeasExt
        NumAddOutputs double = 0 % Number of additional outputs in the system function
        UseSplitModel logical
        minfun = @(X, vecdim) min(X, [], vecdim, 'linear') % min function. The default works for 2019a+.
        failedBackward = 0 % Stage at which the backward phase failed, or 0 otherwise 
        % Properties that control verbosity
        DisplayWarnings logical = true
        DisplayProgressbar logical = true
        % Protected dependent properties
        %   These are defined to avoid property initialization order
        %   dependency along with their dependent counterparts
        VFPenFactorsProt double = [];
        StateNameProt string = [];
        ControlNameProt string = [];
        CostNameProt string = [];
        DisplayProt char {mustBeMember(DisplayProt, {'off', 'warn', 'detailed'})} = 'detailed';
    end
    properties (Dependent)
        N_SV double % total size of the state variables grids
        N_CV double % total size of the control variables grids
        UseExoInput logical
        % Protected dependent properties
        %   These are defined to avoid property initialization order
        %   dependency along with their protected counterparts
        VFPenFactors double;
        StateName
        ControlName
        CostName
        Display 
    end
    properties (SetAccess = private)
        % These properties are user-accessible as read-only
        VF   % Value Function
        Version = "1.6.4";
    end
    properties (Transient, Hidden)
        % Deprecated properties or property names
        VFFactors
        VFInitialization
    end

    %% Methods
    % User-accessible methods
    methods
        % Constructor method
        function obj = DynaProg(StateGrid, StateInitial, StateFinal, ControlGrid, Nstages, varargin)
            %Constructor method
            % Positional arguments:
            %   StateGrid
            %   StateInitial
            %   StateFinal
            %   ControlGrid
            %   Nstages
            % Repeating arguments:
            %   SysName
            %  or
            %   SysNameExt
            %   SysNameInt

            % Handle varargin (Set SysName(s), extract nvps)
            % varargin contains either:
            % - sysname
            % - extsysname, intsysname
            % - sysname, name-value pairs
            % - extsysname, intsysname, name-value pairs
            if isempty(varargin)
                error('DynaProg:wrongNumInputs', ['You must specify at '...
                    'least six positional arguments. Check the '...
                    'syntax guide.'])
            else
                if isa(varargin{1}, 'function_handle') && (length(varargin) == 1 || ~isa(varargin{2}, 'function_handle') )
                    obj.UseSplitModel = false;
                    obj.SysName = varargin{1};
                    NameValuePair = varargin(2:end);
                elseif isa(varargin{1}, 'function_handle') && isa(varargin{2}, 'function_handle')
                    obj.UseSplitModel = true;
                    obj.SysNameExt = varargin{1};
                    obj.SysNameInt = varargin{2};
                    NameValuePair = varargin(3:end);
                else
                    error('DynaProg:wrongSysName', ['You must specify one '...
                        'system function as the sixth positional argument or '...
                        'two system functions as the sixth and seventh '...
                        'arguments. Specify them as function handle(s).'])
                end
            end

            % Set other mandatory arguments
            obj.StateGrid = StateGrid;
            obj.StateInitial = StateInitial;
            obj.StateFinal = StateFinal;
            obj.ControlGrid = ControlGrid;
            obj.Nstages = Nstages;
            % Optional inputs (nvps)
            nvpNames = NameValuePair(1:2:end);
            nvpValues = NameValuePair(2:2:end);
            % Check that the NameValuePair cell array is even
            if mod(length(NameValuePair), 2) ~= 0
                error('DynaProg:oddNVP', ['Additional settings must be '...
                    'specified as name-value pairs.'])
            end
            % Check that the NameValuePair names are string or char array
            for n=1:length(nvpNames)
                if ~ischar(nvpNames{n}) && ~(isstring(nvpNames{n}) && isscalar(nvpNames{n}))
                    error('DynaProg:notStringNVP', ['Additional settings must be '...
                        'specified as name-value pairs. Specify each '...
                        'property name as a string or character array '...
                        'followed by its value.'])
                end
            end
            % Assign the NVP properties
            for n = 1:length(nvpNames)
                obj.(nvpNames{n}) = nvpValues{n};
            end
        end
        % Methods in external files
        obj = run(obj)
        t = plot(obj)
    end

    % Hidden methods (non-user accessible)
    methods (Hidden)
        % Methods in external files
        obj = create_grids(obj)
        obj = create_intVars(obj)
        obj = backward(obj)
        obj = forward(obj)
        [states_next, stageCost, unFeas, addOutputs] = model_wrapper(obj, state, control, exoInput, IntermediateVars)
        [obj, cv_opt_idx, cost_opt] = updateVF(obj, k, states_next, stageCost, unFeas, vecdim_cv)
        [obj, cv_opt, exoInput, intVars_opt] = optimalControl(obj, k, state_next, stageCost, unfeas, vecdim_cv, intVars)
        obj = check_StateFinal(obj)
        obj = checkModelFun(obj, name, mode)
    end

    %% Set/get methods
    methods
        function obj = set.SysNameInt(obj, name)
            obj.SysNameInt = name;
            obj = obj.checkModelFun(name, 'int');

        end
        function obj = set.SysNameExt(obj, name)
            obj.SysNameExt = name;
            obj = obj.checkModelFun(name, 'ext');
        end
        function obj = set.SysName(obj, name)
            obj.SysName = name;
            obj = obj.checkModelFun(name, 'single');
        end
        function obj = set.StateGrid(obj, StateGrid)
            StateGridErr = false;
            if iscell(StateGrid)
                if any(cellfun(@(x) isnumeric(x) || islogical(x), StateGrid))
                    obj.StateGrid  = cellfun(@(x) x(:), StateGrid, 'UniformOutput', false);
                else
                    StateGridErr = true;
                end
            elseif isnumeric(StateGrid) || islogical(StateGrid)
                if isvector(StateGrid)
                    obj.StateGrid = {StateGrid(:)};
                else
                    StateGridErr = true;
                end
            else
                StateGridErr = true;
            end
            if StateGridErr
                error('DynaProg:invalidStateGrid', ['StateGrid must be a '...
                    'vector (if the state is scalar) or a cell array of '...
                    'vectors (if the state is a vector).']);
            end
            % Shift dimensions
            for n = 1:length(obj.StateGrid)
                obj.StateGrid{n} = shiftdim(obj.StateGrid{n}, -(n-1));
            end
        end
        function obj = set.StateInitial(obj, StateInitial)
            if iscell(StateInitial)
                obj.StateInitial = StateInitial;
            else
                obj.StateInitial = num2cell(StateInitial);
            end
            if any(cellfun(@(x) (~isnumeric(x) && ~islogical(x)) || ~isscalar(x), obj.StateInitial))
                error('DynaProg:invalidStateInit', 'StateInitial must be a cell array of scalar values.');
            end
        end
        function obj = set.StateFinal(obj, StateFinal)
            if iscell(StateFinal)
                obj.StateFinal = StateFinal;
            elseif isempty(StateFinal)
                % do nothing
            else
                obj.StateFinal = {StateFinal};
            end
            if any(cellfun(@(x) (~isnumeric(x) && ~islogical(x)) || (length(x)~=2 && ~isempty(x)), obj.StateFinal))
                error('DynaProg:invalidStateFinal', 'StateFinal must be a cell array of two-element numeric vectors or empty values.');
            end
            if any(cellfun(@(x) ~isempty(x) && x(2) < x(1), obj.StateFinal))
                error('DynaProg:wrongOrderStateFinal', 'Each element of StateFinal must contain the lower and upper bound, in this order.');
            end
            if any(cellfun(@(x) ~isempty(x) && (isnan(x(1)) || isnan(x(2))), obj.StateFinal))
                error('DynaProg:nanStateFinal', 'StateFinal does not accept NaNs.');
            end
        end
        function obj = set.ControlGrid(obj, ControlGrid)
            ControlGridErr = false;
            if iscell(ControlGrid)
                if all(cellfun(@(x) isnumeric(x) || islogical(x), ControlGrid))
                    obj.ControlGrid  = ControlGrid;
                else
                    ControlGridErr = true;
                end
            elseif isnumeric(ControlGrid) || islogical(ControlGrid)
                if isvector(ControlGrid)
                    obj.ControlGrid = {ControlGrid(:)};
                else
                    ControlGridErr = true;
                end
            else
                ControlGridErr = true;
            end
            if ControlGridErr
                error('DynaProg:invalidControlGrid', 'ControlGrid must be a vector (if the there is only one control variable) or a cell array of numeric vectors (if there are more).');
            end
        end
        function N_SV = get.N_SV(obj)
            N_SV = cellfun(@(x) length(x), obj.StateGrid);
        end
        function N_CV = get.N_CV(obj)
            N_CV = cellfun(@(x) length(x), obj.ControlGrid);
        end
        function obj = set.ExogenousInput(obj, ExogenousInput)
            if ~iscell(ExogenousInput)
                if isvector(ExogenousInput)
                    obj.ExogenousInput = ExogenousInput(:);
                else
                    error('DynaProg:exogenousInputNotVectors', 'ExogenousInput must be a cell array of vectors, one for each exogenous input variable.')
                end
            else
                for n = 1:length(ExogenousInput)
                    if ~isvector(ExogenousInput{n})
                        error('DynaProg:exogenousInputNotVectors', 'ExogenousInput must be a cell array of vectors, one for each exogenous input variable.')
                    end
                    obj.ExogenousInput(:,n) = ExogenousInput{n}(:);
                end
            end
        end
        function useExoInput = get.UseExoInput(obj)
            if isempty(obj.ExogenousInput)
                useExoInput = false;
            else
                useExoInput = true;
            end
        end
        function VFPenFactors = get.VFPenFactors(obj)
            % If unspecified, define based on the state grid bounds
            if isempty(obj.VFPenFactorsProt)
                % Use a smaller penalty if Level Set is used
                penaltyFactor = ( 1 * obj.UseLevelSet + 10 * ~obj.UseLevelSet );
                for n = 1:length(obj.N_SV)
                    VFPenFactors(n) = abs(obj.StateGrid{n}(end) - obj.StateGrid{n}(1)) * penaltyFactor;
                end
            else
                VFPenFactors = obj.VFPenFactorsProt;
            end
        end
        function obj = set.VFPenFactors(obj, VFPenFactors)
            % Validate size
            if length(VFPenFactors) ~= length(obj.N_SV)
                error('DynaProg:wrongSizeVFPenFactors', 'VFPenFactors must be a numeric array specifying one value for each of the state variables.')
            end
            obj.VFPenFactorsProt = VFPenFactors;
        end
        function stateName = get.StateName(obj)
            % If unspecified (left empty), set StateName to ['x_1', 'x_2', ...]
            if isempty(obj.StateNameProt)
                stateName = string(arrayfun(@(x) ['x_' num2str(x)], 1:length(obj.StateGrid), 'UniformOutput', false));
            else
                stateName = obj.StateNameProt;
            end
        end
        function obj = set.StateName(obj, stateName)
            % If it is a char array, convert it into a cell array
            if ischar(stateName)
                stateName = {stateName};
            end
            % Check the size of StateName
            if length(stateName) ~= length(obj.StateGrid)
                error('DynaProg:wrongSizeStateName', 'StateName must be a string array (or cell array of character vectors) specifying one string for each of the state variables.')
            end
            obj.StateNameProt = stateName;
        end
        function controlName = get.ControlName(obj)
            % If unspecified (left empty), set ControlName to ['u_1', 'u_2', ...]
            if isempty(obj.ControlNameProt)
                controlName = string(arrayfun(@(x) ['u_' num2str(x)], 1:length(obj.ControlGrid), 'UniformOutput', false));
            else
                controlName = obj.ControlNameProt;
            end
        end
        function obj = set.ControlName(obj, controlName)
            % If it is a char array, convert it into a cell array
            if ischar(controlName)
                controlName = {controlName};
            end
            % Check the size of ControlName
            if length(controlName) ~= length(obj.ControlGrid)
                error('DynaProg:wrongSizeControlNames', 'ControlNames must be a string array (or cell array of character vectors) specifying one string for each of the control variables.')
            end
            obj.ControlNameProt = controlName;
        end
        function costName = get.CostName(obj)
            % If unspecified (left empty), set ControlName to 'Cumulative
            % cost'
            if isempty(obj.CostNameProt)
                costName = 'Cumulative cost';
            else
                costName = obj.CostNameProt;
            end
        end
        function obj = set.CostName(obj, costName)
            % If unspecified (left empty), set ControlName to 'Cumulative
            % cost'
            if isempty(costName)
                obj.CostNameProt = 'Cumulative cost';
            else
                obj.CostNameProt = costName;
            end
        end
        function display = get.Display(obj)
            display = obj.DisplayProt;
        end
        function obj = set.Display(obj, Display)
            obj.DisplayProt = Display;
            obj.DisplayProgressbar = strcmp(Display, 'detailed');
            obj.DisplayWarnings = ismember(obj.Display, {'warn', 'detailed'});
        end
        function obj = set.VFInitialization(obj, VFPenalty)
            obj.VFPenalty = VFPenalty;
            warning("Property 'VFInitialization' is deprecated. Use 'VFPenalty' instead.")
        end
        function VFPenalty = get.VFInitialization(obj)
            VFPenalty = obj.VFPenalty;
            warning("Property 'VFInitialization' is deprecated. Use 'VFPenalty' instead.")
        end
        function obj = set.VFFactors(obj, VFPenFactors)
            obj.VFPenFactors = VFPenFactors;
            warning("Property 'VFFactors' is deprecated. Use 'VFPenFactors' instead.")
        end
        function VFPenFactors = get.VFFactors(obj)
            VFPenFactors = obj.VFPenFactors;
            warning("Property 'VFFactors' is deprecated. Use 'VFPenFactors' instead.")
        end
    end
    
    %% Static methods
    methods (Static)
        function [A, linear_index] = min_compatibility(A, vecdim)
            %min_compatibilty min function for releases 2018b and lower
            % Returns the same as
            %   min(A, [], vecdim, 'linear') from 2019a+
            
            ind = cell(1, length(vecdim));
            N_SV = (ndims(A) - length(vecdim));
            subs = cell(1, ndims(A));
            dims = size(A);
            for n = length(ind):-1:1
                [A, ind{n}] = min(A, [], vecdim(n));
            end
            
            if N_SV == 1
                subs{1} = (1:dims(1))';
            elseif N_SV > 1
                [subs{1:N_SV}] = ind2sub(dims, reshape(1:prod(dims(1:N_SV)), dims(1:N_SV)));
            end
            
            subs{N_SV+1} = ind{1};
            for n=1:(length(ind)-1)
                subs{N_SV+n+1} = ind{n+1}( sub2ind( size(ind{n+1}), subs{1:n+N_SV} ) );
            end
            
            linear_index = sub2ind(dims, subs{:});
        end
    end

end