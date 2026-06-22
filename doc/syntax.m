%[text] %[text:anchor:T_305D4575] # DynaProg
%[text] Solve multi-stage deterministic decision problems
%%
%[text] %[text:anchor:H_F09D08F5] ## Syntax
%[text] ```matlabCodeExample
%[text] prob = DynaProg(StateGrid, StateInitial, StateFinal, ControlGrid, NStages, SysName)
%[text] prob = DynaProg(StateGrid, StateInitial, [], ControlGrid, NStages, SysName)
%[text] prob = DynaProg(StateGrid, StateInitial, StateFinal, ControlGrid, NStages, ExtSysName, IntSysName)
%[text] prob = DynaProg(__, Name, Value)
%[text] 
%[text] prob = run(prob)
%[text] plot(prob)
%[text] ```
%%
%[text] %[text:anchor:H_31B35881] ## Description
%[text] [`prob`](internal:H_9FE2499E) `= DynaProg(`[`StateGrid`](internal:H_E1616ABE)`,` [`StateInitial`](internal:H_B28539B7)`,` [`StateFinal`](internal:H_A944F2B8)`,` [`ControlGrid`](internal:H_4CD65B8C)`,` [`NStages`](internal:H_DDBB8435)`,` [`SysName`](internal:H_ED3B20E2)`)` creates the basic problem structure `prob`.
%[text] [`prob`](internal:H_9FE2499E) `= DynaProg(`[`StateGrid`](internal:H_E1616ABE)`,` [`StateInitial`](internal:H_B28539B7)`, [],` [`ControlGrid`](internal:H_4CD65B8C)`,` [`NStages`](internal:H_DDBB8435)`,` [`SysName`](internal:H_ED3B20E2)`)`  creates the basic problem structure `prob`, with no contraints on the final state.
%[text] [`prob`](internal:H_9FE2499E) `= DynaProg(`[`StateGrid`](internal:H_E1616ABE)`,` [`StateInitial`](internal:H_B28539B7)`,` [`StateFinal`](internal:H_A944F2B8)`,` [`ControlGrid`](internal:H_4CD65B8C)`,` [`NStages`](internal:H_DDBB8435)`,` [`ExtSysName`](internal:H_88B8643C)`,` [`IntSysName`](internal:H_038C78B9)`)` creates the basic problem structure prob, using the configuration matrices approach.
%[text] [`prob`](internal:H_9FE2499E) `= DynaProg(__,` [`Name`](internal:H_2F2B9A2D)`,` [`Value`](internal:H_2F2B9A2D)`)` specifies additional properties and information with parameter/value pairs.
%[text] [`prob`](internal:H_9FE2499E) `= run(prob)` runs the optimization algorithm on an existing problem structure and stores the results in the problem structure `prob`.
%[text] `plot(prob)` visualizes results of an existing problem structure containing optimization results.
%%
%[text] %[text:anchor:H_7F17859C] ## Input Arguments
%[text] %[text:anchor:H_E1616ABE] #### `StateGrid`
%[text] `StateGrid` is a cell array of numeric vectors, where each vector defines the discretized grid for a state variable.
%[text] %[text:anchor:H_B28539B7] #### `StateInitial`
%[text] `StateInitial` is a cell array of scalar values, where each value defines the initial value for a state variable.
%[text] %[text:anchor:H_A944F2B8] #### `StateFinal`
%[text] `StateFinal` is a cell array of two-element vectors, where each vector defines lower and upper constraints for the final value for each state variable. Each vector must contain the lower and upper bound, in this order.
%[text] %[text:anchor:H_4CD65B8C] #### `ControlGrid`
%[text] `ControlGrid` is a cell array of numeric vectors, where each vector defines the discretized grid for a control variable.    
%[text] %[text:anchor:H_DDBB8435] #### `NStages` 
%[text] `NStages` is the number of stages of the optimization problem.
%[text] %[text:anchor:H_ED3B20E2] #### `SysName`
%[text] `SysName` is a function handle to the model function. The model function must return the updated state value and the stage cost as a function of the current state value and the control variables. Additionally, it can accept an exogenous input and it can return unfeasibilities.The structure of the model function must be:
%[text]  `[x_next, stageCost, unfeas] = SysName(x, u, w);`
%[text] where:
%[text] - `x` is a cell array, where each cell contains the value for a state variable.
%[text] - `u` is a cell array, where each cell contains the value for a control variable.
%[text] - `w` (optional) is a cell array, where each cell contains the value for an exogenous input variable. If not needed, replace with a tilde (~) in the function signature. See also [Set up a basic problem](file:./set_up_a_basic_problem.mlx).
%[text] - `x_next` is a cell array, where each cell contains the value for the updated state variable.
%[text] - `stageCost` is a numeric variable, containing the stage cost.
%[text] - `unfeas` is either an empty array a boolean variable, set to true for unfeasible (not allowed) combinations of state variables, control variables and exogenous inputs. If not needed, return an empty array (\[\]) in its place. \
%[text] %[text:anchor:H_88B8643C] #### `ExtSysName`
%[text] `ExtSysName` is a function handle to the external model function. The external model function must return the intermediate variables as a function of the control variables. Additionally, it can accept exogenous inputs and it can return unfeasibilities. The structure of the external model function must be:    
%[text]  `[m, unfeas] = ExtSysName(u, w);`
%[text] where:
%[text] `u` is a cell array, where each cell contains the value for a control variable. 
%[text] `w` (optional) is a cell array, where each cell contains the value for an exogenous input variable.
%[text] `m` is a cell array, where each cell contains the value for an intermediate variable.
%[text] `unfeas` (optional) is a boolean variable, set to true for unfeasible (not allowed) combinations of state variables, control variables and exogenous inputs.
%[text] %[text:anchor:H_038C78B9] #### `IntSysName`
%[text] `IntSysName` is a function handle to the internal model function. The internal model function must return the updated state value and the stage cost as a function of the current state value, the control variables and the intermediate variables. Additionally, it can accept an exogenous input and it can return unfeasibilities (both are optional).
%[text] The structure of the model function must be:
%[text] `[x_next, stageCost, unfeas] = IntSysName(x, u, m, w);`
%[text] where:
%[text] - `x` is a cell array, where each cell contains the value for a state variable.
%[text] - `u` is a cell array, where each cell contains the value for a control variable.
%[text] - `m` is a cell array, where each cell contains the value for an intermediate variable.
%[text] - `w` (optional) is a cell array, where each cell contains the value for an exogenous input variable.
%[text] - `x_next` is a cell array, where each cell contains the value for the updated state variable.
%[text] - `stageCost` is a numeric variable, containing the stage cost.
%[text] - `unfeas` (optional) is a boolean variable, set to true for unfeasible (not allowed) combinations of state variables, control variables and exogenous inputs. \
%[text] %[text:anchor:H_2F2B9A2D] ### Name-Value Pair Arguments
%[text] Optional arguments. Specify each argument with a pair of arguments `Name, Value. Name` is one of the argument names specified below and it must be written inside quotes.  `Value` is the value you want to specify for that setting. You can specify any number of name-value pair arguments, in any order.
%[text] %[text:anchor:H_C96247C3] #### **`'ExogenousInput'`**
%[text] Specify exogenous inputs required for your model in a cell array of numeric vectors. Each vector must have the same length as the number of stages of the optimization problem.
%[text] %[text:anchor:H_AA34F1F5] #### **`'UseLevelSet'`**
%[text] Enable Level-Set DP. Defaults to `false` if unspecified.
%[text] #### **`'ForwardMode'`**
%[text] Specify the algorithm to be used in the forward phase. Valid options are `'valueBased'` (default) or `'policyBased'.` See the [Theory](http://./theory.mlx) for more details.
%[text] #### `'SafeMode'`
%[text] Enable Safe Mode (see the [documentation](file:./troubleshooting.mlx)).
%[text] #### `'StoreControlMap'`
%[text] Store the optimal control variables as a function of state for each stage. This cv map is returned in `prob.ControlMap` as a $NC${"altText":"NC"}- by - $N\_{stages}$ cell array, $NC${"altText":"NC"}being the number of control variables. Each cell is a  $N\_{SV\_1} \\times N\_{SV\_2} \\times ... \\times N\_{SV\_\_{NS}} ${"altText":"N\_{SV\_1} x N\_{SV\_2} x ... x N\_{SV\_\_{NS}}"} array which contain the optimal value for the correspondant control variable and stage as a function of the state variables.
%[text] #### **`'StoreValueFunction'`**
%[text] Store the value function. If `true`, running the problem also stores the value function for each stage in a cell array VF.
%[text] %[text:anchor:H_46FC9527] #### `'StateName'`
%[text] Specify state variables names in a string array. 
%[text] %[text:anchor:H_AF88898F] #### `'ControlName'`
%[text] Specify state variables names in a string array.
%[text] %[text:anchor:H_E4B2F10F] #### **`'CostName'`**
%[text] Specify the cumulative cost name as a string.
%[text] %[text:anchor:H_340BB385] #### **`'Time'`**
%[text] Specify time instead of stages. This property is only used in the plots produced with the plot method. It does not affect the optimization.
%[text] #### **`'TerminalCost'`**
%[text] Specify a terminal cost $F(x\_N)$ as a function handle. This property defaults to a null cost. See [Specify a terminal cost](file:./terminal_cost.mlx) and [this example](file:./example_terminal_cost.mlx) for more usage information or the [terminal cost section](file:./theory.mlx:H_6C2F7696) in the theory guide for more information.
%[text] #### **`'VFPenalty'`**
%[text] %[text:anchor:M_B56D103C]  Specify how final state values outside of the final state constraints bounds should be penalized. 
%[text] - Set to `'rift'` to penalize state values outside of the bounds specified as `StateFinal` penalize with a cost equal to infinity. 
%[text] - Set to `'linear'` to penalize them with a cost term proportional to the distance from the bounds.
%[text] - Set to `'none'` to disable any penalization. This may be useful if you want to only use your own terminal cost specified with the `TerminalCost` property. \
%[text] This property defaults to `'linear'` if `'UseLevelSet'` is set to `true` and it is set to `'rift'` otherwise. 
%[text] Read the [terminal cost section](file:./theory.mlx:H_6C2F7696) in the theory guide for more information.
%[text] #### **`'VFPenFactors'`**
%[text] Specify the proportionality factors for the linear terminal cost. This property is only used if VFPenalty is set to `'linear'`. Specify as a numeric array with one factor for each state variable. Read the [terminal cost section](file:./theory.mlx:H_6C2F7696) in the theory guide for more information.
%[text] #### **`'myInf'`**
%[text] Specify penalty term for the terminal cost when  VFPenalty is set to `'rift'`. This term is used to penalize deviations of the terminal state when it is outside of the terminal state bounds. This term should be larger than your expected total cost by at least one or two orders of magnitude. However, it may cause numerical issues if you set it too large.
%[text] #### **`'EnforceStateGrid'`**
%[text] Set a constraint on the state variables so that they do not exceed the state grids. This constraint is needed if your system allows to strongly exceed the state bounduaries over a single stage (simulation step). Defaults to `true` if unspecified.
%[text] #### **`'Display'`**
%[text] Adjust the level of display in the command window. Defaults to `'detailed'` if unspecified.
%[text] - Set to `'off'` to display no output.
%[text] - Set to `'warn'` to displays only warnings and hide the progress bar.
%[text] - Set to `'detailed'` to display both warnings and the progress bar. \
%%
%[text] ## Output Arguments
%[text] %[text:anchor:H_9FE2499E] #### `prob`
%[text] Problem structure, containing all the information required to set up and run the optimization. After `run` is used on a problem structure to run the optmization algorithm, it also contains the optimization results. Most importantly:
%[text] - `prob.StateProfile` is a cell array, with one cell for each state variable. Each cell contains a vector with the values for that state variable at each stage.
%[text] - `prob.ControlProfile` is a cell array, with one cell for each control variable. Each cell contains a vector with the values for that control variable at each stage.
%[text] - `prob.CostProfile` is a vector containing the stageCost at each stage (not including the terminal state cost).
%[text] - `prob.AddOutputsProfile` is a cell array, with one cell for each additional output. Each cell contains a vector with the values for that additional output at each stage.
%[text] - `prob.totalCost` is the total evaluated cost. If the optimization failed, `prob.totalCost` is set to `inf`. \
%%
%[text]{"align":"right"} *Contact: federico.miretti@polito.it*

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
