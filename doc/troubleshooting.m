%[text] # Troubleshooting your model and Safe Mode
%[text] In order to obtain the optimal solution to your control problem, DynaProg uses your system and cost function(s) with n-dimensional arrays. Some operations may produce unintended results. In order to avoid these issues, you can either follow this [troubleshooting guide](internal:H_D4C99F2E) or [enable Safe Mode](internal:H_F4B9C99C).
%[text] %[text:anchor:H_D4C99F2E] ## Troubleshooting your model function
%[text] In the backward phase of the control optimization algorithm, DynaProg uses your system and cost function to evaluate updated state values, stage costs and unfeasibilities. When it does so, it passes  state and control variables passes as an input in the form of n-dimensional arrays. In these arrays, each dimension is assigned to one specific state or control variable. 
%[text] Consider an optimization problem with two state variables and two control variables, which are discretized on computational grids having $N\_{x\_1}${"altText":"N\_(x\_1)"}, $N\_{x\_2}${"altText":"N\_(x\_2)"}, $N\_{u\_1}${"altText":"N\_(u\_1)"} and $N\_{u\_2}${"altText":"N\_(u\_2)"} elements respectively. The system and cost function is defined as:
function [x_next, stage_cost, unfeas] = sysfun(x, u, ~)
    
a = u{1} .* u{2};
x_next{1} = x{1} + a .* dt;
x_next{2} = x{1} + (x{1} - x{2}) .* dt;

stage_cost = u{1};

unfeas = [];
end
%[text] In the backward phase, DynaProg passes as inputs:
%[text] - `x{1}` as an array of size $N\_{x\_1} \\times 1 \\times 1 \\times 1 ${"altText":"N\_(x\_1) x 1 x 1 x 1"},
%[text] - `x{2}` as an array of size $1 \\times N\_{x\_2} \\times 1 \\times 1 ${"altText":"1 x N\_(x\_2) x 1 x 1"},
%[text] - `u{1}` as an array of size $1 \\times 1 \\times N\_{u\_1} \\times 1 ${"altText":"1 x 1 x N\_(u\_1) x 1"},
%[text] - `u{1}` as an array of size $1 \\times 1 \\times 1 \\times N\_{u\_2}${"altText":"1 x 1 x 1 x N\_(u\_2)"}. \
%[text] As MATLAB evaluates each line of the system and cost function, it will expand the size of the resulting array accordingly. Therefore:
%[text] - `a` will be an array of size $1 \\times 1 \\times N\_{u\_1} \\times N\_{u\_2} ${"altText":"1 x 1 x N\_(u\_1) x N\_(u\_2)"}, as it is only dependent on the control variables,
%[text] - `x_next{1}` will be an array of size $N\_{x\_1} \\times 1 \\times N\_{u\_1} \\times N\_{u\_2} ${"altText":"N\_(x\_1) x 1 x N\_(u\_1) x N\_(u\_2)"}, as it is only dependent on $x\_1${"altText":"x\_1"}, $u\_1${"altText":"u\_1"} and $u\_2${"altText":"u\_2"},
%[text] - `x_next{2}` will be an array of size $N\_{x\_1} \\times N\_{x\_2} \\times N\_{u\_1} \\times N\_{u\_2} ${"altText":"N\_(x\_1) x N\_(x\_2) x N\_(u\_1) x N\_(u\_2)"}, as it is only dependent on all state and control variables,
%[text] - `stageCost` will be an array of size $1 \\times 1 \\times N\_{u\_1} \\times 1${"altText":"1 x 1 x N\_(u\_1) x 1"}, as it is only dependent on $u\_1${"altText":"u\_1"}. \
%[text] For this reason, you must ensure that all operations that you use in your system and cost function are consistent in their output when operating on n-dimensional arrays.
%[text] Example of such operations are:
%[text] - Element-wise sum (`+`), subtraction (`-`), product (`.*`), division (`./`) and power (`.^`).
%[text] - Indexing with an nd-array with at least two non-singleton dimensions.
%[text] - Logical indexing if the indexed array and the index array have the same size.
%[text] - Functions that operate elementwise and preserve the input size. \
%[text] Example of operations that are inconsistent are:
%[text] - Matrix operations such as matrix multiplication (`*`), "division" (`\` and `/`) and power (`^`).
%[text] - Indexing with an nd-array with at only one one non-singleton dimensions (i.e. a vector).
%[text] - Logical indexing if the indexed array and the index array do not have the same size. \
%[text] To see some examples of troublesome operations, and how to fix them, see [this example](file:./example_troubleshooting.mlx).
%[text] %[text:anchor:H_F4B9C99C] ## Safe Mode
%[text] Safe Mode is an alternative way of dealing with the system and cost function which is meant to automatically handle potentially troublesome operations, at the cost of increased optimization time.
%[text] By default, Safe Mode is disabled. To enable Safe Mode, specify it in the constructor after the mandatory positional arguments:
prob = DynaProg(StateGrid, StateInitial, StateFinal, ControlGrid, Nstages, SysName, 'SafeMode', true);
%[text] or set the `SafeMode` property to `true` on an existing problem structure.
prob.SafeMode = true;
%[text] Safe Mode forces DynaProg to always pass the state variables, control variables, exogenous inputs and intermediate variables as full n-dimensional arrays having size $N\_{x\_1} \\times N\_{x\_2} \\times  ... \\times  N\_{x\_{NX}} \\times N\_{u\_1} \\times N\_{u\_2} \\times  ... \\times N\_{u\_{NU}}$  where $N\_{x\_i}${"altText":"N\_x\_i"} and $N\_{u\_j}${"altText":"N\_u\_j"} are the number of grid points that discretize the state variable $x\_i${"altText":"x\_i"} and the control variable $u\_j${"altText":"u\_j"} and$NX${"altText":"NX"} and $NU${"altText":"NU"}are the number of state and control variables. 
%[text] This effectively avoids all errors related to array indexing, but it means that many repetitive, unnecessary computations are performed. Note, however, that Safe Mode may not prevent all sorts of errors that may arise when you run the model function.

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
