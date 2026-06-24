%[text] # Specify a terminal cost
%[text] To specify a terminal cost, you need to create a function handle and use that to set the `TerminalCost` property.
%[text] ### Basic usage
%[text] The handle must take one input, which represents a cell array containing the state variables $x\_N$.
prob.TerminalCost = @(x) myTerminalCost(x);
%[text] Note that you can pass any additional parameters to your function handle. For example, if your terminal cost must somehow depend on a desired terminal state `x_final_constraint`, which you defined somewhere in your script, you can define `TerminalCost` as
x_final = 10;
prob.TerminalCost = @(x) myTerminalCost(x, x_final);
%[text] ### Terminal cost and terminal state constraint
%[text] Note that, if a terminal state is specified, DynaProg adds by default a penalty term $\\Psi(x\_N)$ to the terminal cost in order to enforce the corresponding terminal state constraints. In other words, the cost-to-go at stage $N\n$ is initialized to:
%[text] $V\_N(x\_N) = F(x\_N) + \\Psi(x\_N)$.
%[text] Where $F(x\_N)$ is the cost specified with the `TerminalCost` property (defaults to empty). If you want to disable this feature to gain full control on your terminal cost, you must set the `VFPenalty` propery to `'none'`: 
prob.VFPenalty = 'none';
%[text] which will set $\\Psi(x\_N)$ to zero.
%[text] ### Example
%[text] See the [double-integrator with custom terminal cost example](file:./example_terminal_cost.mlx).

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
