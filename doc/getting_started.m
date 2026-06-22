%[text] %[text:anchor:T_240B1D76] # Getting started
%[text] Solve finite horizon multi-stage deterministic decision problems
%%
%[text] DynaProg provides a flexible tool to solve a finite horizon multi-stage deterministic decision problem, which is a problem where a decision must be made at each stage for a system that evolves through a finite number of stages, minimizing the total cost incurred.
%[text] This class of problems is described by a dynamic system and a cost function.
%[text] The dynamic system usually takes the form:
%[text] $x\_{k+1} = f(x\_k, u\_k) \\qquad  k = 0,1,...,N-1,\n${"altText":"x\_(k+1) = f(x\_k, u\_k), for k = 0,1,...,N-1,"}
%[text] where 
%[text] - $k${"altText":"k"} indexes the stages (which may represent discrete time intervals or something more abstract),
%[text] - $x\_k\n${"altText":"x\_k"} characterizes the state of the system at stage $k${"altText":"k"},
%[text] - $u\_k${"altText":"u\_k"} characterizes the control variables (the decision) to be taken at stage $k${"altText":"k"},
%[text] - $N${"altText":"N"} is the total number of stages. \
%[text] The cost function takes the form:
%[text] $J(x\_0,u\_0,...,u\_{N-1}) = g\_N(x\_N) + \\sum\_{k=0}^{N-1} g(x\_k, u\_k),${"altText":"J(x\_0,u\_0,...,u\_(N-1)) = g\_N(x\_N) + SUM(from k=0 to N-1) of g(x\_k, u\_k),"}
%[text] where 
%[text] - $g(x\_k,u\_k)${"altText":"g(x\_k,u\_k)"} is the stage cost, that is the cost incurred for each stage, and
%[text] - $g\_N(x\_N)${"altText":"g\_N(x\_N)"} is some terminal cost. \
%[text] Obtaining the optimal solution of this optimization problem is therefore to obtain the minimum total cost  
%[text] $J^\*(x\_0) = \\min\_{u\_k, k=0,\\\> ...,\\\> N-1} J(x\_0,u\_0,...,u\_{N-1})${"altText":"J^\*(x\_0) = MIN over (u\_k, k=0, ..., N-1) of J(x\_0,u\_0,...,u\_{N-1})"}
%[text] and the optimal control sequence $u^\*\_0,\\\>...,\\\>u^\*\_{N-1}${"altText":"u\*\_0, ..., u\*\_(N-1)"} that minimizes it.
%[text] With DynaProg, you can model the dynamic system and a cost function as a MATLAB function in an m-file, and easily obtain the optimal solution, without having to develop your own implementation of a Dynamic Programming algorithm.
%[text] Additionally, you can specify constraints on the state and control variables and/or add an exogenous input which can influence the system dynamics.
%%
%[text] ## [Set up a basic problem](file:./set_up_a_basic_problem.mlx)
%[text] Set up and solve a basic problem with DynaProg.
%[text] ## [Specify contraints on the state and control variables](file:./specify_constraints.mlx)
%[text] Use unfeasibilities in your system to add constraints on the state and control variables.
%[text] ## [Specify a terminal cost](file:./terminal_cost.mlx)
%[text] Use a custom function to define a terminal cost.
%[text] ## [Handle exogenous inputs](file:./handle_exogenous_inputs.mlx)
%[text] Consider an exogenous input influencing the system dynamics.
%[text] ## [Use additional inputs](file:./use_additional_inputs.mlx)
%[text] Use additional inputs to your model to avoid repeating evaluations.
%[text] ## [Use additional outputs](file:./use_additional_outputs.mlx)
%[text] Return additional simulation outputs.
%[text] ## [Splitting the model function](file:./splitting_the_model.mlx)
%[text] Speed up the solution for complex models by splitting the system and cost function.
%[text] ## [Troubleshooting your model and Safe Mode](file:./troubleshooting.mlx)
%[text] Speed up the solution for complex models by splitting the system and cost function.
%%
%[text]{"align":"right"} *Contact: federico.miretti@polito.it*

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":35}
%---
