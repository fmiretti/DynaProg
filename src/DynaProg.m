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
    %       myInf: if VFInitialization is set to 'rift', myInf defines the
    %           penalty cost for unfeasible terminal states.
    %       EnforceStateGrid: set a constraint so that the state variables
    %           do not exceed the state grids. Defaults to true.
    %     # Terminal Cost (Value Function initialization) settings
    %       VFInitialization: specify how final state values outside of 
    %           the final state constraints bounds should be penalized. Set
    %           to 'rift' to penalize them with a myInf term. Set to
    %           'linear' to penalize them with a term proportional to the 
    %           distance from the bounds. Set to 'manual' to use a custom 
    %           terminal cost and define it with TerminalCost. Note that in
    %           this case DynaProg will not attempt to enforce STATEFINAL.
    %       TerminalCost: define a custom terminal cost. Its size must be
    %           composed by the lenghts of the state variable grids (i.e.
    %           length(x1_grid) x length(x2_grid) x ... x length(xn_grid))
    %       VFFactors: if VFInitialization is set to 'linear', VFFactors 
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
        VFInitialization char {mustBeMember(VFInitialization, {'rift', 'linear', 'auto', 'manual'})} = 'auto';
        LevelSetInitialization char = [];
        TerminalCost double = [];
        % Results 
        StateProfile
        ControlProfile
        CostProfile
        AddOutputsProfile
        ControlMap
    end
    properties (Access = protected)
        % Computational grids, value and level-set function
        StateCombGrid
        StateGridVect % State grids as N_SV_n-by-1 vectors, needed to create VF interpolants
        ControlCombGrid
        ControlFullGrid % Full CV grid
        ControlFullShiftedGrid
        LevelSet % Level-Set function
        IntermediateVars
        unFeasExt
        NumAddOutputs double = 0 % Number of additional outputs in the system function
        UseSplitModel logical
        minfun = @(X, vecdim) min(X, [], vecdim, 'linear') % min function. The default works for 2019a+.
        % Protected dependent properties
        %   These are defined to avoid property initialization order
        %   dependency along with their dependent counterparts
        VFFactorsProt double = [];
        StateNameProt string = [];
        ControlNameProt string = [];
        CostNameProt string = [];
    end
    properties (Dependent)
        N_SV double % total size of the state variables grids
        N_CV double % total size of the control variables grids
        UseExoInput logical
        % Protected dependent properties
        %   These are defined to avoid property initialization order
        %   dependency along with their protected counterparts
        VFFactors double;
        StateName
        ControlName
        CostName
    end
    properties (SetAccess = private)
        VF   % Value Function
        Version = "1.6";
    end
    
    methods
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
        
        function obj = run(obj)
            % Check the MATLAB version
            if verLessThan('matlab','9.6') % aka 2019a
                obj.minfun =  @(X, vecdim) obj.min_compatibility(X, vecdim);
            end
            % Run the optimization algorithm
            obj = create_grids(obj);
            if obj.UseSplitModel
                obj = create_intVars(obj);
            end
            obj = backward(obj);
            obj = forward(obj);
            % Remove VF if its storage was not required
            if ~obj.StoreValueFunction
                obj.VF = [];
            end
            % Reset warnings
            warning('on', 'DynaProg:failedCostToGoUpdate')
        end
        
        function obj = create_grids(obj)
            % Create computational grids, initialize terminal value
            % function and level-set function
            
            % checks
            if ~isempty(obj.StateGrid)
                obj = check_StateFinal(obj);
                if length(obj.StateFinal) ~= length(obj.StateGrid)
                    error('DynaProg:wrongSizeStateFinal', ['You must specify '...
                        'one final condition for each of the state variables.']);
                end
            end
            if length(obj.StateInitial) ~= length(obj.StateGrid)
                error('DynaProg:wrongSizeStateInit', ['You must specify one '...
                    'initial condition for each of the state variables.']);
            end
            % Set the VF initialization method if unspecified
            if strcmp(obj.VFInitialization, 'auto')
                if isempty(obj.TerminalCost)
                    if obj.UseLevelSet
                        obj.VFInitialization = 'linear';
                    else
                        obj.VFInitialization = 'rift';
                    end
                else
                    obj.VFInitialization = 'manual';
                end
            end
             % Checks on the VF initialization method
            if ~isempty(obj.TerminalCost) && ~strcmp(obj.VFInitialization, 'manual')
                warning('DynaProg:ignoredTerminalCost', 'You specified a terminal cost with TerminalCost but VFInitialization is not set to ''manual''. Ignoring your terminal cost.' )
            end
            if isempty(obj.TerminalCost) && strcmp(obj.VFInitialization, 'manual')
                warning('DynaProg:emptyTerminalCost', 'You set VFInitialization to ''manual'' but you did not specify a terminal cost. Set VFInitialization to a valid string or specify a terminal cost with TerminalCost.' )
            end
             % Set the Level Set initialization method if unspecified
            if isempty(obj.LevelSetInitialization)
                obj.LevelSetInitialization = 'linear';
            end
            
            % State grids as N_SV_n-by-1 vectors, needed to create VF interpolants
            obj.StateGridVect = cellfun(@(x) x(:), obj.StateGrid, 'UniformOutput', false);
            
            % Combined SV grids
            obj.StateCombGrid = cell(1, length(obj.N_SV));
            obj.ControlCombGrid = cell(1, length(obj.N_CV));
            [obj.StateCombGrid{:}, obj.ControlCombGrid{:}] = ndgrid(obj.StateGrid{:}, obj.ControlGrid{:});
                
            % Full CV grid
            obj.ControlFullGrid = cell(1, length(obj.ControlGrid));
            [obj.ControlFullGrid{:}] = ndgrid(obj.ControlGrid{:});
            for n = 1:length(obj.N_CV)
                obj.ControlGrid{n} = obj.ControlGrid{n}(:);
                obj.ControlGrid{n} = shiftdim(obj.ControlGrid{n}, -length(obj.StateGrid) - (n-1));
                obj.ControlFullShiftedGrid{n} = shiftdim(obj.ControlFullGrid{n}, -length(obj.StateGrid));
            end
            % Initialize terminal VF
            StateFullGrid = cell(1, length(obj.StateGrid));
            [StateFullGrid{:}] = ndgrid(obj.StateGrid{:});
            VFN = zeros(size(StateFullGrid{1}));
            if ~isempty(obj.StateFinal)
                switch obj.VFInitialization
                    case 'linear' %VFN is proportional to the distance from the target set
                        for n = 1:length(obj.StateGrid)
                            if ~isempty(obj.StateFinal{n})
                                VFN = VFN + obj.VFFactors(n) .* ( max(obj.StateFinal{n}(1)-StateFullGrid{n}, 0) + max(StateFullGrid{n}-obj.StateFinal{n}(2), 0) );
                            end
                        end
                        VFN(isinf(VFN)) = 0;
                    case 'rift' % VFN is inf outside the target set
                        for n=1:length(obj.StateGrid)
                            if ~isempty(obj.StateFinal{n})
                                VFN( StateFullGrid{n} > obj.StateFinal{n}(2) ...
                                    | StateFullGrid{n} < obj.StateFinal{n}(1)) = obj.myInf;
                            end
                        end
                    case 'manual'
                        VFN = obj.TerminalCost;
                        % Check user-supplied terminal cost
                        if ~isequal(size(VFN), obj.N_SV)
                            error('DynaProg:wrongSizeTerminalCost', strjoin({'The terminal cost you provided has wrong size. It should be', sprintf('%dx', obj.N_SV), sprintf('\b\b (the lengths of the state grids).'), '\n'}))
                        end
                end
            end
            % Allocate cell array for the VF interpolants and create the interpolant
            obj.VF = cell(1, obj.Nstages+1);
            obj.VF{end} = griddedInterpolant(obj.StateGridVect, VFN, ...
                'linear');
            % Initialize Level Set function
            if obj.UseLevelSet
                LevelSetN = cell(1,length(obj.StateGrid));
                for n=1:length(obj.StateGrid)
                    if ~isempty(obj.StateFinal{n})
                        LevelSetN{n} = max( obj.StateFinal{n}(1) - StateFullGrid{n}, ...
                            StateFullGrid{n} - obj.StateFinal{n}(2) );
                    else
                        LevelSetN{n} = zeros(size(StateFullGrid{n}));
                    end
                end
                LevelSetN = cat(length(LevelSetN)+1,LevelSetN{:});
                LevelSetN = max(LevelSetN,[],ndims(LevelSetN));
                
                % Allocate cell array for the LevelSet interpolants and create the interpolant
                obj.LevelSet = cell(1, obj.Nstages+1);
                obj.LevelSet{end} = griddedInterpolant(obj.StateGridVect, LevelSetN, ...
                    'linear');
                
            end
        end
        
        function obj = create_intVars(obj)
            % Run external model to create intermediate variables
            exoInput = cell(1, size(obj.ExogenousInput, 2));
            if obj.SafeMode
                control = obj.ControlFullGrid;
            else
                control = obj.ControlGrid;
            end
            for k = 1:obj.Nstages
                % Create exogenous inputs
                if obj.UseExoInput
                    currentExoInput = obj.ExogenousInput(k,:);
                    for n = 1:length(currentExoInput)
                        if obj.SafeMode
                            exoInput{n} = currentExoInput(n).*ones(size(obj.ControlFullGrid{1}));
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
        
        function obj = backward(obj)
            % Run the optimization algorithm backward phase 
            if obj.SafeMode
                state = obj.StateCombGrid;
                control = obj.ControlCombGrid;
            else
                state = obj.StateGrid;
                control = obj.ControlGrid;
            end
                
            % Progress Bar
            fprintf('DP backward progress:    ')
            % Vector dimensions corresponding to CVs
            vecdim_cv = (length(obj.N_SV)+1):(length(obj.N_CV)+length(obj.N_SV));
            
            % Backward Loop
            for k = obj.Nstages:-1:1
                % Progress Bar
                fprintf('%s%2d %%', ones(1,4)*8, floor((obj.Nstages-k)/obj.Nstages*100));
                if obj.UseSplitModel
                    if obj.SafeMode
                        % Expand the intermediate variables to the combined grid
                        intVars = cellfun(@(x) repmat(shiftdim(x, -length(obj.N_SV)), [obj.N_SV ones(size(obj.N_CV))]), obj.IntermediateVars{k}, 'UniformOutput', false);
                        unfeasExt = repmat(shiftdim(obj.unFeasExt{k}, -length(obj.N_SV)), [obj.N_SV ones(size(obj.N_CV))]);
                    else
                        intVars = obj.IntermediateVars{k};
                        unfeasExt = obj.unFeasExt{k};
                    end
                else
                    intVars = [];
                end
                % Create exogenous inputs
                if obj.UseExoInput
                    currentExoInput = obj.ExogenousInput(k,:);
                    for n = 1:length(currentExoInput)
                        if obj.SafeMode
                            exoInput{n} = currentExoInput(n).*ones(size(obj.StateCombGrid{1}));
                        else
                            exoInput{n} = currentExoInput(n);
                        end
                    end
                else
                    exoInput = [];
                end
                
                % State update
                [states_next, stageCost, unFeasInt] = model_wrapper(obj, state, control, exoInput, intVars);
                unFeasInt = logical(unFeasInt);
                
                % Model output checks
                if k == obj.Nstages
                    if iscell(stageCost)
                        error('DynaProg:wrongFormatStageCost', ['The stage '...
                            'cost must be returned as a numeric type, not a cell.'])
                    elseif ~isnumeric(stageCost)
                        error('DynaProg:wrongFormatStageCost', ['The stage '...
                            'cost must be returned as a numeric type.'])
                    end
                end
                
                % feasibility: include external model unfeasibility
                if obj.UseSplitModel
                    unFeas = unFeasInt | unfeasExt;
                else
                    unFeas = unFeasInt;
                end
                if ~obj.SafeMode
                    % Expand updated states and unfeas to the combined grid
                    for n = 1:length(states_next)
                        states_next{n} = states_next{n} + zeros([obj.N_SV obj.N_CV]);
                    end
                    stageCost = stageCost + zeros([obj.N_SV obj.N_CV]);
                    unFeas = unFeas | false([obj.N_SV obj.N_CV]);
                end
                % Enforce state grids
                if obj.EnforceStateGrid
                    for n = 1:length(obj.N_SV)
                        unFeas(states_next{n} > obj.StateGrid{n}(end) | states_next{n} < obj.StateGrid{n}(1)) = obj.myInf;
                    end
                end
                
                % Update Level-Set function
                if obj.UseLevelSet
                    % Read L(k+1)
                    LevelSet_next =  obj.LevelSet{k+1}(states_next{:});
                    % set LevelSet_next to inf for the unfeasible CVs
                    LevelSet_next(unFeas) = obj.myInf;
                    % Update level-set function and find L-minimizing CVs
                    [LevelSetValue, MinLevelSetCV] = obj.minfun(LevelSet_next, vecdim_cv);
                    % Check if the set of reachable CVs U^R(x_k) (CVs that
                    % lead to a reachable state) is empty. Note: U^R(x_k) is
                    % a subset of U(x_k) (set of feasible CVs).
                    isempty_UR = LevelSetValue>0;
                    % Level set failure
                    if all(LevelSetValue>0)
                        %                         keyboard;
                    end
                    % Construct L approximating function for the current timestep
                    obj.LevelSet{k} = griddedInterpolant(obj.StateGridVect, ...
                        LevelSetValue, 'linear');
                end
                
                % Read VF(k+1)
                VF_next =  obj.VF{k+1}(states_next{:});
                cost = stageCost + VF_next;
                % Set cost-to-go to inf for the unfeasible/unreachable CVs
                cost(unFeas) = obj.myInf;
                % Find optimal control as a function of the current state
                [cost_opt, cv_opt] = obj.minfun(cost, vecdim_cv);
                
                if obj.UseLevelSet
                    % For those state grid points where no feasible cv was found
                    % (isempty_UR), calculate the VF based on the cv that minimizes
                    % the level-set function.
                    cost_MinLevelSetCV = cost(MinLevelSetCV);
                    cost_opt(isempty_UR) = cost_MinLevelSetCV(isempty_UR);
                    if obj.StoreControlMap
                        cv_opt(isempty_UR) = MinLevelSetCV(isempty_UR);
                    end
                end
                
                % Construct VF approximating function for the current timestep
                obj.VF{k} = griddedInterpolant(obj.StateGridVect, cost_opt, ...
                    'linear');
                
                % Store cv map
                if obj.StoreControlMap
                    for n = 1:length(obj.N_CV)
                        cv_opt_sub = cell(1, ndims(cost));
                        [cv_opt_sub{:}] = ind2sub(size(cost), cv_opt);
                        obj.ControlMap{n,k} = squeeze(obj.ControlGrid{n}(cv_opt_sub{n+length(obj.N_SV)}));
                    end
                end
                
                % Warn the user if there are no feasible trajectories for
                % the tail subproblem
                if all(cost_opt >= obj.myInf) 
                    warning('DynaProg:failedCostToGoUpdate', 'The Cost-to-Go update has failed. Your problem might be overconstrained or the state variables grid might be too coarse.')
                    warning('off', 'DynaProg:failedCostToGoUpdate')
                end
                
            end
            if obj.VF{1}(obj.StateInitial) >= obj.myInf
                fprintf('\n')
                error('DynaProg:failedCostToGoUpdate', 'The Cost-to-Go update has failed: no feasible solution was found. Your problem might be overconstrained or the state variables grid might be too coarse.\n')
            end

            % Progress Bar
            fprintf('%s%2d %%\n', ones(1,4)*8, 100);
        end
        
        function obj = forward(obj)
            % Run the optimization algorithm forward phase
            
            if obj.SafeMode
                control = obj.ControlFullGrid;
                % Vector dimensions corresponding to CVs
                vecdim_cv = 1:length(obj.N_CV);
            else
                control = obj.ControlGrid;
                % Vector dimensions corresponding to CVs
                vecdim_cv = (1:length(obj.N_CV)) + length(obj.N_SV);
            end

            % Progress Bar
            progressbar = true;
            fprintf('DP forward progress:    ')
            % Preallocate profiles
            StateProfileMat = zeros(length(obj.StateGrid), obj.Nstages+1);
            StateProfileMat(:,1) = [obj.StateInitial{:}];
            ControlProfileMat = zeros(length(obj.ControlGrid), obj.Nstages);
            % Initialize the state
            state = obj.StateInitial;
            obj.CostProfile(1) = 0;
            state_next = cell(1, length(obj.N_SV));
            exoInput = cell(1, size(obj.ExogenousInput, 2));
            % Initialize warnings
            unfeasFwdWarn = false;
            unfeasFwdWarnStages = [];

            for k = 1:obj.Nstages
                % Progress Bar
                if progressbar
                    fprintf('%s%2d %%', ones(1,4)*8, floor((k-1)/obj.Nstages*100));
                end

                % Expand current state to the full cv grid
                if obj.SafeMode
                    for n = 1:length(state)
                        state_next{n} = state{n} + zeros(size(obj.ControlFullGrid{1}));
                    end
                else
                    state_next = state;
                end
                % Create intermediate vaiables
                if obj.UseSplitModel
                    intVars = obj.IntermediateVars{k};
                    unfeasExt = obj.unFeasExt{k};
                else
                    intVars = [];
                end
                % Create exogenous inputs
                if obj.UseExoInput
                    exoInput_temp = obj.ExogenousInput(k,:);
                    for n = 1:length(exoInput_temp)
                        if obj.SafeMode
                            exoInput{n} = exoInput_temp(n).*ones(size(obj.ControlFullGrid{1}));
                        else
                            exoInput{n} = exoInput_temp(n);
                        end
                    end
                else
                    exoInput = [];
                end
                
                % Evaluate state update and stage cost
                [state_next, stageCost, unFeas] = model_wrapper(obj, state_next, control, exoInput, intVars);
                unFeas = logical(unFeas);
                if obj.UseSplitModel
                    unFeas = unFeas | unfeasExt;
                end
                
                % Expand updated states and unfeas to the full cv grid
                if ~obj.SafeMode
                    for n = 1:length(state_next)
                        state_next{n} = state_next{n} + zeros([ones(1, length(obj.N_SV)) obj.N_CV]);
                    end
                    stageCost = stageCost + zeros([ones(1, length(obj.N_SV)) obj.N_CV]);
                    unFeas = unFeas | false([ones(1, length(obj.N_SV)) obj.N_CV]);
                end
                
                % Get Level Set-minimizing CV
                if obj.UseLevelSet
                    % Read L(k+1)
                    LevelSet_next = obj.LevelSet{k+1}(state_next{:});
                    % set LevelSet_next to inf for the unfeasible CVs
                    LevelSet_next(unFeas) = obj.myInf;
                    % Determine if U^R(x_k) is empty
                    isempty_UR = all(LevelSet_next(:) > 0);
                    % Find L-minimizing u
                    [~, MinLevelSetCV] = obj.minfun(LevelSet_next, vecdim_cv);
                end
                
                % Read VF(k+1)
                VF_next =  obj.VF{k+1}(state_next{:});
                cost = stageCost + VF_next;
                if obj.UseLevelSet
                    cost(unFeas) = obj.myInf;
                    cost(LevelSet_next > 0) = obj.myInf;
                end
                % Set cost-to-go to inf for the unfeasible/unreachable CVs
                cost(unFeas) = obj.myInf;
                
                % Find optimal control as a function of the current state
                [~, index_opt] = obj.minfun(cost, vecdim_cv);
                if obj.UseLevelSet
                    % If no reachable cv was found (isempty_UR), select the cv that
                    % minimizes the level-set function.
                    index_opt(isempty_UR) = MinLevelSetCV;
                end
                cv_opt =  cellfun(@(x) x(index_opt), obj.ControlFullGrid, 'UniformOutput', false);
                
                % Extract the intermediate variables for the optimal cv
                if obj.UseSplitModel
                    intVars = cellfun(@(x) x .* ones(size(cost)), intVars, 'UniformOutput', false);
                    intVars_opt = cellfun(@(x) x(index_opt), intVars, 'UniformOutput', false);
                else
                    intVars_opt = [];
                end
                % Advance the simulation
                [state, stageCost, unfeas, addout] = model_wrapper(obj, state, cv_opt, exoInput, intVars_opt);
                
                % Check solution validity
                if unfeas
                    unfeasFwdWarnStages(end+1) = k;
                    unfeasFwdWarn = true;
                end

                % Update the profiles
                StateProfileMat(:,k+1) = [state{:}]';
                ControlProfileMat(:,k) = [cv_opt{:}]';
                if ~isempty(addout)
                    for n = 1:length(addout)
                        obj.AddOutputsProfile{n}(k) = addout{n};
                    end
                end
                obj.CostProfile(k+1) = stageCost;
            end
            % Check terminal state constraints
            for n = 1:length(obj.StateFinal)
                if ~isempty(obj.StateFinal{n})
                    if state{n} < obj.StateFinal{n}(1) || state{n} > obj.StateFinal{n}(2)
                        fprintf('\n')
                        warning('DynaProg:failedTerminalState', ['The solution violates your terminal state constraints. Your problem might be overconstrained or the state variables grid might be too coarse.\n' ...
                            'You can try refining the grids, widening the final state constraint bounds or using the Level Set option.\n'])
                        progressbar = false;
                        break
                    end
                end
            end
            % Print information about constraints violation in the fwd run
            if unfeasFwdWarn
                unfeasFwdWarnStages(end+1) = k;
                unfeasFwdWarnStages = unfeasFwdWarnStages(1:min(10, numel(unfeasFwdWarnStages)));
                str = "The solution violates your constraints at stages:\n";
                for n = 1:numel(unfeasFwdWarnStages)
                    str = str + "%d ";
                end
                if numel(unfeasFwdWarnStages) >= 10
                    str = str + "and others more";
                end
                str = str + ".\nYour problem might be overconstrained or the state variables grid might be too coarse.";
                unfeasFwdWarnStages = num2cell(unfeasFwdWarnStages);
                warning('DynaProg:failedForward', str, unfeasFwdWarnStages{:})
            end
            % Store state and control profiles
            obj.StateProfile = num2cell(StateProfileMat,2);
            obj.ControlProfile = num2cell(ControlProfileMat,2);
            % Progress Bar
            if progressbar
                fprintf('%s%2d %%\n', ones(1,4)*8, 100);
            end
        end
        
        function [states_next, stageCost, unFeas, addOutputs] = model_wrapper(obj, state, control, exoInput, IntermediateVars)
            getAddOuts = nargout > 3;
            if getAddOuts
                addOutputs = cell(1, obj.NumAddOutputs);
            end
            % Run the function model
            if obj.UseSplitModel
                if ~getAddOuts && ~obj.UseExoInput
                    [states_next, stageCost, unFeas] = obj.SysNameInt(state, control, IntermediateVars);
                    addOutputs = {};
                elseif  ~getAddOuts && obj.UseExoInput
                    [states_next, stageCost, unFeas] = obj.SysNameInt(state, control, exoInput, IntermediateVars);
                    addOutputs = {};
                elseif getAddOuts && ~obj.UseExoInput
                    [states_next, stageCost, unFeas, addOutputs{:}] = obj.SysNameInt(state, control, IntermediateVars);
                else
                    [states_next, stageCost, unFeas, addOutputs{:}] = obj.SysNameInt(state, control, exoInput, IntermediateVars);
                end
            else
                if ~getAddOuts && ~obj.UseExoInput
                    [states_next, stageCost, unFeas] = obj.SysName(state, control);
                    addOutputs = {};
                elseif  ~getAddOuts && obj.UseExoInput
                    [states_next, stageCost, unFeas] = obj.SysName(state, control, exoInput);
                    addOutputs = {};
                elseif  getAddOuts && ~obj.UseExoInput
                    [states_next, stageCost, unFeas, addOutputs{:}] = obj.SysName(state, control);
                else
                    [states_next, stageCost, unFeas, addOutputs{:}] = obj.SysName(state, control, exoInput);
                end
            end
            if isempty(unFeas)
                unFeas = zeros(size(state{1}));
            end
        end
        
        function obj = check_StateFinal(obj)
            for n = 1:length(obj.StateGrid)
                if ~isempty(obj.StateFinal{n})
                    if obj.UseLevelSet & (obj.StateFinal{n}(1) == obj.StateFinal{n}(2))
                        warning('DynaProg:levelSetFail', 'When the level-set method is enabled, you should not specify the target state as a single point. Make sure you specify distinct lower and upper final state bounds.')
                    else
                        count = sum(obj.StateGrid{n} >= obj.StateFinal{n}(1) & obj.StateGrid{n} <= obj.StateFinal{n}(2));
                        if count == 0
                            warning('DynaProg:tooCoarseStateGrid', ['The state variable grid for SV #' num2str(n) ' is too coarse. At least one of the grid points should lie inside the final state bounds.'])
                        elseif count < 3
                            warning('DynaProg:coarseStateGrid', ['The state variable grid for SV #' num2str(n) ' is very coarse. Consider defining the grid so that at least three of the grid points lie inside the final state bounds, or loosen the final state bounds.'])
                        end
                    end
                end
            end
        end
        
        function t = plot(obj)
            %DynaProg.plot Plot optimization results
            %   plot(prob) plots state variables, control variables and
            %   cumulative cost profiles.
            
            % Check the MATLAB version
            if verLessThan('matlab','9.7') % aka 2019b
                compatibility_mode = true;
            else
                compatibility_mode = false;
            end
            
            if isempty(obj.StateProfile)
                error('DynaProg:notSolved', 'This problem structure does not contain any solution. Run the optimization to produce results.')
            end
            % Set up layout
            ncols = lcm(length(obj.N_SV), length(obj.N_CV));
            sv_span = ncols/length(obj.N_SV);
            cv_span = ncols/length(obj.N_CV);
            
            if compatibility_mode
                subplot(3, ncols, 1);
            else
                t = tiledlayout(3, ncols);
            end
            
            ax = [];
            % Plot SV profiles
            for n = 1:length(obj.N_SV)
                if compatibility_mode
                    ax(end+1) = subplot(3, ncols, [1 + sv_span*(n-1) sv_span*n]);
                else
                    ax(end+1) = nexttile([1 sv_span]); %#ok<*AGROW>
                end
                if isempty(obj.Time)
                    plot(obj.StateProfile{n}, 'LineWidth', 1.5)
                else
                    plot(obj.Time, obj.StateProfile{n}, 'LineWidth', 1.5)
                end
                title(obj.StateName(n))
                axis tight
            end
            % Plot CV profiles
            for n = 1:length(obj.N_CV)
                if compatibility_mode
                    ax(end+1) = subplot(3, ncols, [ncols + 1 + cv_span*(n-1), ncols + cv_span*n]);
                else
                    ax(end+1) = nexttile([1 cv_span]); %#ok<*AGROW>
                end
                if isempty(obj.Time)
                    plot(obj.ControlProfile{n}, 'LineWidth', 1.5)
                else
                    plot(obj.Time(1:end-1), obj.ControlProfile{n}, 'LineWidth', 1.5)
                end
                title(obj.ControlName(n))
                axis tight
            end
            % Plot cumulative cost profile
            if compatibility_mode
                ax(end+1) = subplot(3, ncols, [ncols*2 + 1, ncols*3]);
            else
                ax(end+1) = nexttile([1 ncols]); %#ok<*AGROW>
            end
            if isempty(obj.Time)
                cumCost = cumsum(obj.CostProfile);
                plot(cumCost, 'LineWidth', 1.5)
                xlabel('Stage number')
            else
                cumCost = cumsum(obj.CostProfile);
                plot(obj.Time, cumCost, 'LineWidth', 1.5)
                xlabel('Time')
            end
            title(obj.CostName)
            axis tight
            
            % Finalize            
            if compatibility_mode
                arrayfun(@(x) set(x, 'FontSize', 10), ax)
                t = gcf;
            else
                arrayfun(@(x) set(x, 'FontSize', 10), t.Children)
                t.TileSpacing = 'compact';
                t.Padding = 'compact';
            end
            linkaxes(ax,'x')
        end
        
    end

    % Set/get methods
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
        function VFFactors = get.VFFactors(obj)
            % If unspecified, define based on the state grid bounds
            if isempty(obj.VFFactorsProt)
                for n = 1:length(obj.N_SV)
                    VFFactors(n) = abs(obj.StateGrid{n}(end) - obj.StateGrid{n}(1)) * 10;
                end
            else
                VFFactors = obj.VFFactorsProt;
            end
        end
        function obj = set.VFFactors(obj, VFFactors)
            % Validate size
            if length(VFFactors) ~= length(obj.N_SV)
                error('DynaProg:wrongSizeVFFactors', 'VFFactors must be a numeric array specifying one value for each of the state variables.')
            end
            obj.VFFactorsProt = VFFactors;
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
    end
    
    methods (Access=protected)
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

    end

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