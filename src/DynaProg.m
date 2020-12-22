classdef (CaseInsensitiveProperties=true) DynaProg
    %DynaProg  Create DynaProg problem structure.
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
    %       UseLevelSet: Enable Level-Set DP. Defaults to false if
    %       unspecified.
    %       ExogenousInput: specify exogenous inputs required for your
    %       model in a cell array of numeric vectors. Each vector must have
    %       the same length as the number of stages of the optimization
    %       problem.
    %       SafeMode: enables Safe Mode (see the documentation)
    %       StateName: specify state variables names in a string
    %       array. 
    %       ControlName: specify state variables names in a string
    %       array.
    %       CostName: specify the cumulative cost name as a string.
    %       Time: specify time instead of stages. This property is only
    %       used in the plots produced with the plot method. It does not 
    %       affect the optimization.
    %   
    %   Author: Federico Miretti
    %
    % '<a href="matlab:web html/index.html">Apri qui</a>'
    
    properties
        % Constructor: Main properties
        SysName function_handle
        SysNameExt function_handle
        SysNameInt function_handle
        StateGrid cell
        StateInitial cell
        StateFinal cell = [];
        ControlGrid cell
        Nstages
        % Constructor: Name-Value pairs
        StateName string
        ControlName string
        CostName string
        ExogenousInput
        UseLevelSet logical = false;
        VFInitialization char
        LevelSetInitialization char
        SafeMode = false;
        Time double
        % Results 
        StateProfile
        ControlProfile
        CostProfile
        AddOutputsProfile
    end
    properties (Access = private)
        % Computational grids, value and level-set function
        N_SV double % total size of the state variables grids
        N_CV double % total size of the control variables grids
        StateCombGrid
        ControlCombGrid
        ControlFullGrid % Full CV grid
        ControlFullShiftedGrid
        myInf = 1e6;
        VF   % Value Function
        LevelSet % Level-Set function
        IntermediateVars
        unFeasExt
        UseExoInput
        NumAddOutputs
        UseSplitModel logical
    end
    
    methods
        function obj = DynaProg(StateGrid, StateInitial, StateFinal, ControlGrid, Nstages, SysName, NameValuePair)
            arguments
                StateGrid
                StateInitial
                StateFinal
                ControlGrid
                Nstages
            end
            arguments(Repeating)
                SysName function_handle
            end
            arguments
                NameValuePair.StateName = [];
                NameValuePair.ControlName = [];
                NameValuePair.CostName = [];
                NameValuePair.ExogenousInput = [];
                NameValuePair.UseLevelSet = false;
                NameValuePair.SafeMode = false;
                NameValuePair.VFInitialization = [];
                NameValuePair.LevelSetInitialization = [];
                NameValuePair.Time = [];
            end
            % Set SysName(s)
            switch length(SysName)
                case 1
                    obj.UseSplitModel = false;
                    obj.SysName = SysName{1};
                case 2
                    obj.UseSplitModel = true;
                    obj.SysNameExt = SysName{1};
                    obj.SysNameInt = SysName{2};
                otherwise
                    error("Do not specify more than two model functions.")
            end
            % Set other mandatory arguments
            obj.StateGrid = StateGrid;
            obj.StateInitial = StateInitial;
            obj.StateFinal = StateFinal;
            obj.ControlGrid = ControlGrid;
            obj.Nstages = Nstages;
            % Optional inputs
            obj.StateName = NameValuePair.StateName;
            obj.ControlName = NameValuePair.ControlName;
            obj.CostName = NameValuePair.CostName;
            obj.ExogenousInput = NameValuePair.ExogenousInput;
            obj.UseLevelSet = NameValuePair.UseLevelSet;
            obj.SafeMode = NameValuePair.SafeMode;
            obj.VFInitialization  = NameValuePair.VFInitialization;
            obj.LevelSetInitialization  = NameValuePair.LevelSetInitialization;
            obj.Time = NameValuePair.Time;
            % Set others
            obj.N_SV = cellfun(@(x) length(x), obj.StateGrid);
            obj.N_CV = cellfun(@(x) length(x), obj.ControlGrid);
        end
        
        function obj = run(obj)
            % Run the optimization algorithm
            obj = create_grids(obj);
            if obj.UseSplitModel
                obj = create_intVars(obj);
            end
            obj = backward(obj);
            obj = forward(obj);
        end
        
        function obj = create_grids(obj)
            % Create computational grids, initialize terminal value
            % function and level-set function
            
            % checks
            if ~isempty(obj.StateGrid)
                obj = check_StateFinal(obj);
                if length(obj.StateFinal) ~= length(obj.StateGrid)
                    error('DynaProg:wrongSizeStateFinal', "You must specify one final condition for each of the state variables.");
                end
            end
            if length(obj.StateInitial) ~= length(obj.StateGrid)
                error('DynaProg:wrongSizeStateInit', "You must specify one initial condition for each of the state variables.");
            end
            if isempty(obj.VFInitialization)
                if obj.UseLevelSet
                    obj.VFInitialization = 'linear';
                else
                    obj.VFInitialization = 'ridge';
                end
            end
            if isempty(obj.LevelSetInitialization)
                obj.LevelSetInitialization = 'linear';
            end
            
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
                                VFN = VFN + max(obj.StateFinal{n}(1)-StateFullGrid{n}, 0) + max(StateFullGrid{n}-obj.StateFinal{n}(2), 0);
                            end
                        end
                        VFN(isinf(VFN)) = 0;
                        if ~obj.UseLevelSet
                            VFN = VFN.*1e1;
                        end
                    case 'ridge' % VFN is inf outside the target set
                        for n=1:length(obj.StateGrid)
                            if ~isempty(obj.StateFinal{n})
                                VFN( StateFullGrid{n} > obj.StateFinal{n}(2) ...
                                    | StateFullGrid{n} < obj.StateFinal{n}(1)) = obj.myInf;
                            end
                        end
                end
            end
            % Allocate cell array for the VF interpolants and create the interpolant
            obj.VF = cell(1, obj.Nstages+1);
            obj.VF{end} = griddedInterpolant(obj.StateGrid, VFN, ...
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
                LevelSetN = cat(ndims(LevelSetN)+1,LevelSetN{:});
                LevelSetN = max(LevelSetN,[],ndims(LevelSetN));
                
                % Allocate cell array for the LevelSet interpolants and create the interpolant
                obj.LevelSet = cell(1, obj.Nstages+1);
                obj.LevelSet{end} = griddedInterpolant(obj.StateGrid, LevelSetN, ...
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
                sizes = cellfun(@(x) size(x), obj.IntermediateVars{k}, "UniformOutput", false);
                if obj.SafeMode
                    wrong_size = ~cellfun(@(x) isequal(x, obj.N_CV), sizes);
                else
                    for n = 1:length(sizes)
                        sizes_sv = sizes{n}(1:length(obj.N_SV));
                        sizes_cv = sizes{n}(length(obj.N_SV)+1:end);
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
            fprintf("DP backward progress:    ")
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
                
                % Update Level-Set function
                if obj.UseLevelSet
                    % Read L(k+1)
                    LevelSet_next =  obj.LevelSet{k+1}(states_next{:});
                    % set LevelSet_next to inf for the unfeasible CVs
                    LevelSet_next(unFeas) = obj.myInf;
                    % Update level-set function and find L-minimizing CVs
                    [LevelSetValue, MinLevelSetCV] = min(LevelSet_next, [], vecdim_cv, 'linear');
                    % Check if the set of reachable CVs U^R(x_k) (CVs that
                    % lead to a reachable state) is empty. Note: U^R(x_k) is
                    % a subset of U(x_k) (set of feasible CVs).
                    isempty_UR = LevelSetValue>0;
                    % Level set failure
                    if all(LevelSetValue>0)
                        %                         keyboard;
                    end
                    % Construct L approximating function for the current timestep
                    obj.LevelSet{k} = griddedInterpolant(obj.StateGrid, ...
                        LevelSetValue, 'linear');
                end
                
                % Read VF(k+1)
                VF_next =  obj.VF{k+1}(states_next{:});
                cost = stageCost + VF_next;
                % Set cost-to-go to inf for the unfeasible/unreachable CVs
                cost(unFeas) = obj.myInf;
                % Find optimal control as a function of the current state
                [cost_opt, ~] = min(cost, [], vecdim_cv, 'linear');
                if obj.UseLevelSet
                    % For those state grid points where no feasible cv was found
                    % (isempty_UR), calculate the VF based on the cv that minimizes
                    % the level-set function.
                    cost_MinLevelSetCV = cost(MinLevelSetCV);
                    cost_opt(isempty_UR) = cost_MinLevelSetCV(isempty_UR);
                end
                
                % Construct VF approximating function for the current timestep
                obj.VF{k} = griddedInterpolant(obj.StateGrid, cost_opt, ...
                    'linear');
                % Break the dp bwd run if there are no feasible trajectories for
                % the tail subproblem
                if all(cost_opt >= obj.myInf)
                    warning('The Cost-to-Go update has failed. Your problem might be overconstrained or the state variables grid might be too coarse.')
                end
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
            fprintf("DP forward progress:    ")
            % Preallocate profiles
            StateProfileMat = zeros(length(obj.StateGrid), obj.Nstages+1);
            StateProfileMat(:,1) = [obj.StateInitial{:}];
            ControlProfileMat = zeros(length(obj.ControlGrid), obj.Nstages);
            % Initialize the state
            state = obj.StateInitial;
            obj.CostProfile(1) = 0;
            state_next = cell(1, length(obj.N_SV));
            exoInput = cell(1, size(obj.ExogenousInput, 2));
            for k = 1:obj.Nstages
                % Progress Bar
                fprintf('%s%2d %%', ones(1,4)*8, floor((k-1)/obj.Nstages*100));
                % Expand current state
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
                % Get Level Set-minimizing CV
                if obj.UseLevelSet
                    % Read L(k+1)
                    LevelSet_next = obj.LevelSet{k+1}(state_next{:});
                    % set LevelSet_next to inf for the unfeasible CVs
                    LevelSet_next(unFeas) = obj.myInf;
                    % Determine if U^R(x_k) is empty
                    isempty_UR = all(LevelSet_next > 0, 'all');
                    % Find L-minimizing u
                    [~, MinLevelSetCV] = min(LevelSet_next, [], 'all', 'linear');
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
                [~, index_opt] = min(cost, [], vecdim_cv, 'linear');
                if obj.UseLevelSet
                    % If no reachable cv was found (isempty_UR), select the cv that
                    % minimizes the level-set function.
                    index_opt(isempty_UR) = MinLevelSetCV;
                end
                cv_opt =  cellfun(@(x) x(index_opt), obj.ControlFullGrid, 'UniformOutput', false);
                if obj.UseSplitModel
                    intVars_opt = cellfun(@(x) x(index_opt), intVars, 'UniformOutput', false);
                else
                    intVars_opt = [];
                end
                % Advance the simulation
                [state, stageCost, unFeas, addout] = model_wrapper(obj, state, cv_opt, exoInput, intVars_opt);
                StateProfileMat(:,k+1) = [state{:}]';
                ControlProfileMat(:,k) = [cv_opt{:}]';
                if ~isempty(addout)
                    for n = 1:length(addout)
                        obj.AddOutputsProfile{n}(k) = addout{n};
                    end
                end
                % optional
                obj.CostProfile(k+1) = stageCost;
            end
            obj.StateProfile = num2cell(StateProfileMat,2);
            obj.ControlProfile = num2cell(ControlProfileMat,2);
            % Progress Bar
            fprintf('%s%2d %%\n', ones(1,4)*8, 100);
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
            
            if isempty(obj.StateProfile)
                error('DynaProg:notSolved', 'This problem structure does not contain any solution. Run the optimization to produce results.')
            end
            % Set up layout
            ncols = lcm(length(obj.N_SV), length(obj.N_CV));
            sv_span = ncols/length(obj.N_SV);
            cv_span = ncols/length(obj.N_CV);
            t = tiledlayout(3, ncols);
            ax = [];
            % Plot SV profiles
            for n = 1:length(obj.N_SV)
                ax(end+1) = nexttile([1 sv_span]); %#ok<*AGROW>
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
                ax(end+1) = nexttile([1 cv_span]);
                if isempty(obj.Time)
                    plot(obj.ControlProfile{n}, 'LineWidth', 1.5)
                else
                    plot(obj.Time(1:end-1), obj.ControlProfile{n}, 'LineWidth', 1.5)
                end
                title(obj.ControlName(n))
                axis tight
            end
            % Plot cumulative cost profile
            ax(end+1) = nexttile([1 ncols]);
            if isempty(obj.Time)
                cumCost = cumsum(obj.CostProfile);
                plot(cumCost, 'LineWidth', 1.5)
            else
                cumCost = cumtrapz(obj.Time, obj.CostProfile);
                plot(obj.Time, cumCost, 'LineWidth', 1.5)
            end
            title(obj.CostName)
            axis tight
            % Finalize
            if isempty(obj.Time)
                xlabel(t, 'Stage number')
            else
                xlabel(t, 'Time')
            end
            arrayfun(@(x) set(x, 'FontSize', 10), t.Children)
            t.TileSpacing = 'compact';
            t.Padding = 'compact';
            linkaxes(ax,'x')
        end
        
    end
    
    % Set methods
    methods
        function obj = set.SysNameInt(obj, name)
            try
                obj.SysNameInt = name;
            catch ME
                modelNameErr(ME)
            end
            obj = obj.checkAddOutputs(name);
        end
        function obj = set.SysNameExt(obj, name)
            try
                obj.SysNameExt = name;
            catch ME
                modelNameErr(ME)
            end
        end
        function obj = set.SysName(obj, name)
            try
                obj.SysName = name;
            catch ME
                modelNameErr(ME)
            end
            obj = obj.checkAddOutputs(name);
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
                error('DynaProg:invalidStateGrid', "StateGrid must be a vector (if the state is scalar) or a cell array of vectors (if the state is a vector).");
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
                error('DynaProg:invalidStateInit', "StateInitial must be a cell array of scalar values.");
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
                error('DynaProg:invalidStateFinal', "StateFinal must be a cell array of two-element numeric vectors or empty values.");
            end
            if any(cellfun(@(x) ~isempty(x) && x(2) < x(1), obj.StateFinal))
                error('DynaProg:wrongOrderStateFinal', "Each element of StateFinal must contain the lower and upper bound, in this order.");
            end
            if any(cellfun(@(x) ~isempty(x) && (isnan(x(1)) || isnan(x(2))), obj.StateFinal))
                error('DynaProg:nanStateFinal', "StateFinal does not accept NaNs.");
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
                error('DynaProg:invalidControlGrid', "ControlGrid must be a vector (if the there is only one control variable) or a cell array of numeric vectors (if there are more).");
            end
        end
        function obj = set.StateName(obj, StateName)
            if isempty(StateName)
                obj.StateName = string(arrayfun(@(x) ['x_' num2str(x)], 1:length(obj.StateGrid), 'UniformOutput', false));
            else
                obj.StateName = StateName;
            end
            if length(obj.StateName) ~= length(obj.StateGrid)
                error('DynaProg:wrongSizeStateName', 'StateName must be a string array specifying one string for each of the state variables.')
            end
        end
        function obj = set.ControlName(obj, ControlName)
            if isempty(ControlName)
                obj.ControlName = string(arrayfun(@(x) ['u_' num2str(x)], 1:length(obj.ControlGrid), 'UniformOutput', false));
            else
                obj.ControlName = ControlName;
            end
            if length(obj.ControlName) ~= length(obj.ControlGrid)
                error('DynaProg:wrongSizeControlNames', 'ControlNames must be a string array specifying one string for each of the control variables.')
            end
        end
        function obj = set.CostName(obj, CostName)
            if isempty(CostName)
                obj.CostName = "Cumulative cost";
            else
                obj.CostName = CostName;
            end
        end
        function obj = set.ExogenousInput(obj, ExogenousInput)
            if isempty(ExogenousInput)
                obj.UseExoInput = false;
            elseif ~iscell(ExogenousInput)
                if isvector(ExogenousInput)
                    obj.ExogenousInput = ExogenousInput(:);
                    obj.UseExoInput = true;
                else
                    error('DynaProg:exogenousInputNotVectors', 'ExogenousInput must be a cell array of vectors, one for each exogenous input variable.')
                end
            else
                for n = 1:length(ExogenousInput)
                    if ~isvector(ExogenousInput{n})
                        error('DynaProg:exogenousInputNotVectors', 'ExogenousInput must be a cell array of vectors, one for each exogenous input variable.')
                    end
                    obj.ExogenousInput(:,n) = ExogenousInput{n}(:);
                    obj.UseExoInput = true;
                end
            end
        end
    end
    
    methods (Access = private)
        function obj = checkAddOutputs(obj, name)
            info = functions(name);
            if strcmp(info.type, 'anonymous')
                func_name = regexp(info.function, '\)\w*\(', 'match');
                func_name = func_name{1}(2:end-1);
            else
                func_name = name;
            end
            if nargout(func_name) < 3
                error('DynaProg:invalidModelOutput', "The function model has a wrong number of outputs.\n")
            end
            obj.NumAddOutputs = nargout(func_name) - 3;
        end
    end
    
    methods (Static)
        function modelNameErr(ME)
            switch ME.identifier
                case 'MATLAB:narginout:functionDoesnotExist'
                    error(['Model function ' char(name) ' not found. Check the filename and path.'])
                otherwise
                    rethrow(ME)
            end
        end
    end
end