%[text] # Handle exogenous inputs
%[text] The system dynamics and the stage cost fundamentally depend on the current state of the system itself and the control variables. 
%[text] - system dynamics    $x\_{k+1} = f(x\_k, u\_k) \\qquad  k = 0,1,...,N-1,\n${"altText":"x\_(k+1) = f(x\_k, u\_k), for k = 0,1,...,N-1,"}
%[text] - stage cost               $g(x\_k,u\_k)${"altText":"g(x\_k,u\_k)"}  \
%[text] However, it is also possible to introduce a dependence of the state dynamics and stage cost on another signal wich is known a priori for each stage of the decision problem. In the context of DynaProg, this type of signal is called an exogenous input ($w\_k${"altText":"w\_k"}). Thus, the sytem dynamics and stage cost can take the form
%[text] - system dynamics    $x\_{k+1} = f(x\_k, u\_k, w\_k) \\qquad  k = 0,1,...,N-1,\n${"altText":"x\_(k+1) = f(x\_k, u\_k, w\_k), for k = 0,1,...,N-1,"}
%[text] - stage cost               $g(x\_k,u\_k,w\_k)${"altText":"g(x\_k, u\_k, w\_k)"}.  \
%[text] #### Include an exogenous input the system and cost function
%[text] Create the system and cost function with signature:
[x_next, stage_cost, unfeas] = sysfun(x, u, w)
%[text] where `w` is used to pass exogenous input to your function. Similarly to `x` and `u`, it must be treated as a cell array, where each cell can contain a separate exogenous input.
%[text] #### Include the exogenous input in the problem structure
%[text] In the script where you set up and solve the problem, create the DynaProg problem by also including the Name-Value pair argument `'ExogenousInput'`, followed by a cell array of vectors which contains the value of your exogenous inputs for each stage of the problem.
w{1} = exo_input_1;
w{2} = exo_input_2;
prob = DynaProg(___, 'ExogenousInput', w);
%[text] Each vector contained in `w` should have the same length as the number of stages of the problem.
%[text] 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
