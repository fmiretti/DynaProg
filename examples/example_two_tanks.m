%Two-Tank system example
% Taken from: P. Elbert, S. Ebbesen, and L. Guzzella, “Implementation of
% Dynamic Programming for n-Dimensional Optimal Control Problems With Final
% State Constraints,” IEEE Transactions on Control Systems Technology, vol.
% 21, no. 3, pp. 924–931, May 2013, doi: 10.1109/TCST.2012.2190935.

%% Set up the problem
clear
% State variables grids
SVnames = ["Tank #1 level", "Tank #2 level"];
x1_grid = 0:0.01:1;
x2_grid = 0:0.01:1;
x_grid = {x1_grid, x2_grid};
% Initial state 
x1_init = 0;
x2_init = 0;
x_init = {x1_init, x2_init};
% Final state constraints
x1_final = [0.5 inf];
x2_final = [0.5 inf];
x_final = {x1_final, x2_final};
% Control variables grids
CVnames = ["Throttle position", "Direction valve position"];
u1_grid = 0:0.05:1;
u2_grid = 0:0.05:1;
u_grid = {u1_grid, u2_grid};
% Number of stages (time intervals)
Nint = 200;
time = 0:0.01:2;

% Create DynaProg object
prob = DynaProg(x_grid, x_init, x_final, u_grid, Nint, @two_tanks);

%% Solve and visualize results
prob.UseLevelSet = true;
prob = run(prob);

prob.plot;