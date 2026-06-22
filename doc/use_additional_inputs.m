%[text] # Use additional inputs
%[text] Using additional inputs is strongly recommended for all those situations where the system and cost function requires to access large amounts of data which does not change from one stage to the other (such as in the [HEV example](file:./example_hev.mlx)).
%[text] The reason for this is that, while the system and cost function is very flexbile in allowing to define a custom optimization problem, it should be noted that the optimization algorithm evaluates the system and cost function a very high number of times (see [Notes on computational time](file:./theory.mlx)).
%[text] To use additional inputs in your problem, you must:
%[text] - Specify any additional input in your function signature
%[text] - Parametrize the function to be used in constructing the problem structure. \
%[text] ## Include additional inputs in the system and cost function
%[text] Create the system and cost function with signature:
[x_next, stage_cost, unfeas] = sysfun(x, u, w, addInput1, .., addInputN)
%[text] You can specify any number of additional inputs, but they always must be specified starting from the fourth input as the first three inputs are reserved to the state variables, control variables and exogenous inputs. 
%[text] If you wish to specify additional inputs but you do not wish to specify exogenous inputs, replace the third input argument with a tilde `~`:
[x_next, stage_cost, unfeas] = sysfun(x, u, ~, addInput1, .., addInputN)
%[text] ## Parametrize the system and cost function
%[text] In the script where you set up and solve the problem, create or load the data to be passed as additional inputs. 
%[text] Then, when creating the problem structure, specify SysName as a function handle (or anonymous function) that parametrizes your system and cost function with your additional inputs. For example, if you are also using exogenous inputs, create your anonymous function as:
SysName = @(x, u, w) sysfun(x, u, w, addInput1, .., addInputN)
%[text] or, if you are not using exogenous inputs:
SysName = @(x, u, w) sysfun(x, u, [], addInput1, .., addInputN)
%[text] ## Differences between exogenous inputs and additional inputs
%[text] Both the exogenous and additional inputs are important ways to reduce unnecessary computations in the system and cost function, which is critical to the overall optimization time. Both can be used to define variables that influence the state dynamics but are not influenced by the state and control variables themselves. However, there is one notable difference between them:
%[text] - Exogenous inputs are useful to define quantities whose value is stage-dependent. For example, it can easily represent external, time-dependent signals.
%[text] - Additional inputs are useful to define quantities which are stage-indepentent. For example, it may be used to define constant physical characteristics of the system. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
