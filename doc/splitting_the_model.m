%[text] # Splitting the model function
%[text] Splitting the system and cost function may significantly reduce the optimization time for complex models where many operations are evaluated involving several state and control varialbles.
%[text] In particular, it is effective if at least part of the operations performed by the system and cost function do not involve the state variables.
%[text] Consider the standard form for the state dynamics and stage cost in DynaProg:
%[text] - system dynamics    $x\_{k+1} = f(x\_k, u\_k, w\_k) \\qquad  k = 0,1,...,N-1,\n${"altText":"x\_(k+1) = f(x\_k, u\_k, w\_k), for k = 0,1,...,N-1,"}
%[text] - stage cost               $g(x\_k,u\_k,w\_k)${"altText":"g(x\_k, u\_k, w\_k)"}.  \
%[text] Splitting the model function means to split the system dynamics and the stage cost in two distinct phases.
%[text] - system dynamics    $x\_{k+1} = f\_{int}(x\_k,u\_k, w\_k,  f\_{ext}(u\_k, w\_k) )\\qquad  k = 0,1,...,N-1,\n${"altText":"x\_(k+1) = f(x\_k, u\_k, w\_k), for k = 0,1,...,N-1,"}
%[text] - stage cost               $g\_{int}(x\_k,u\_k,w\_k, g\_{ext}(u\_k,w\_k))${"altText":"g(x\_k, u\_k, w\_k)"}.  \
%[text] The *external* model performs all calculations in which the values of the state variables are not involved. The *internal* model performs the remaining calculations. 
%[text] The reason why the external model is named so is because it can be run before the backward phase (see [theory](file:./theory.mlx)), on a different (much smaller in size) computational grid which does not involve the state variables. Performing part of the calculations on this reduced computational grid is the advantage of this approach.
%[text] In practice, you will:
%[text] - Create an external function to generate some *intermediate* *variables*,
%[text] - create an internal function to evaluate the updated system state and the stage cost based on the state variables, control variables,  exogenous inputs (optionally) and the intermediate variables. \
%[text] ## Create the external system and cost function
%[text] Create the external system and cost function with signature:
function [intVar, unfeas] = sysfun_ext(u, w)
   ...
end
%[text] Your function should accept as the first input the control variables `u.` The second input can be used to pass an [exogenous input](file:./handle_exogenous_inputs.mlx) to the system function. If your system model does not make use of exogenous inputs, you can replace the second input in your function signature with a tilde (~).
%[text] Also, the function should return two outputs: the first output contains the intermediate variables. The intermediate variables are all relevant quantities that you will need to use in the internal function to ultimately evaluate the state update and the stage cost. Treat the intermediate variables as a cell array, with each cell containing one intermediate variable.
%[text] The second output is used to define unfeasibility. Note that the unfeasibilties that you specify in the external function will be combined with those specified in the internal model function, so that both will  be taken into account. If you do not wish to dpecify unfeasibilty, you must still specify a second output in your function signature, but this should return an empty array: you should include anywhere in your function the line `unfeas = [];.`
%[text] You can specify any number of [additional inputs](file:./use_additional_inputs.mlx), but you should not specify any [additional output](file:./use_additional_outputs.mlx) as those must be included in the internal function.
function [intVar, unfeas] = sysfun_ext(u, w, addInput1, .., addInputN)
   ...
end
%[text] ## Create the internal system and cost function
[x_new, stageCost, unfeas] = sysfun_ext(x, u, w, intVar)
%[text] Your function should accept four inputs: the first three are the state variables, the control variables and the exogenous inputs`.` The fourth input will contain the intermediate variables that your external function returns; remember to treat it as a cell array. Once again, the exogenous inputs can be replaced by a tilde `~` if not needed.
%[text] As for the outputs, your internal function should follow the same rules as for regular [system and cost functions](file:./set_up_a_basic_problem.mlx). The first three outputs are the updated state variables, the stage cost and unfeasibility. 
%[text] You can also specify any number of [additional inputs](file:./use_additional_inputs.mlx) and [additional outputs](file:./use_additional_outputs.mlx).
%[text] ## Set up and solve the optimization problem
%[text] Set up the state and and control variables grids, initial conditions and terminal constraints [as usual](file:./set_up_a_basic_problem.mlx).
%[text] Pass a function handle to the external function as the sixth input argument and a function handle to the internal function as the seventh argument.
prob = DynaProg(StateGrid, StateInitial, StateFinal, ControlGrid, Nstages, ...
 @sysfun_ext, @sysfun_int);
%[text] You can specify any value-pair argument after the two function handles.
prob = DynaProg(StateGrid, StateInitial, StateFinal, ControlGrid, Nstages, ...
 @sysfun_ext, @sysfun_int, 'ExogenousInput', w, 'StateName', SVnames);
%[text] To pass [additional inputs](file:./use_additional_inputs.mlx) to your external and internal functions, parametrize them as you normally would, making sure you specify the correct inputs for the function handles. 
prob = DynaProg(StateGrid, StateInitial, StateFinal, ControlGrid, Nstages, ...
 @(u, w) sysfun_ext(u, w, addInput1, .., addInputN), @(x, u, w, m) sysfun_int(x, u, w, m, addInput1, .., addInputN), 'ExogenousInput', w);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
