%HEV example #2 using Safe Mode
% Control strategy optimization of a p2 parallel HEV passenger car with
% three-way catalyst
% Data source: ADVISOR software (2003, "ADVISOR". National Renewable Energy
% Laboratory. Available: https://sourceforge.net/projects/adv-vehicle-sim/,
% accessed on: Nov 2020)

%% Set up the problem
clear
% State variables grids
SVnames = {'SOC', 'TWC Temperature'};
x1_grid = 0.4:0.001:0.7;
x2_grid = 200:10:600;
x_grid = {x1_grid, x2_grid};
% Initial state 
x1_init = 0.6;
x2_init = 300;
x_init = {x1_init, x2_init};
% Final state constraints
x1_final = [0.599 0.601];
x2_final = [];
x_final = {x1_final, x2_final};
% Control variables grids
CVnames = {'Gear Number', 'Torque split'};
u1_grid = [1 2 3 4 5];
u2_grid = -1:0.1:1;
u_grid = {u1_grid, u2_grid};

% Load a drive cycle
load UDDS % contains velocity and time vectors
dt = time_s(2) - time_s(1);
% Create exogenous input
w{1} = speed_kmh./3.6; 
w{2} = [diff(w{1})/dt; 0]; 

% Number of stages (time intervals)
Nint = length(time_s);

% Generate and store vehicle data
[veh, fd, gb, eng, em, batt, twc] = hev_data();

% Create DynaProg object (using Safe Mode)
prob = DynaProg(x_grid, x_init, x_final, u_grid, Nint, ...
 @(u, w) hev_ext_naive(u, w, veh, fd, gb, eng, em), @(x, u, w, m) hev_int_naive(x, u, w, m, batt, twc), ...
 'ExogenousInput', w, 'SafeMode', true);

%% Solve and visualize results
% Solve the problem
prob = run(prob);

% Add time vector to use in the plot
prob.Time = [time_s; time_s(end)+dt];
% Add SV, CV and cost names to be used in the plot
prob.StateName = SVnames;
prob.ControlName = CVnames;

% Plot results
fc_grams = prob.AddOutputsProfile{1};
hcTailpipeFlwRate = [prob.AddOutputsProfile{2}.hc];
coTailpipeFlwRate = [prob.AddOutputsProfile{2}.co];
noxTailpipeFlwRate = [prob.AddOutputsProfile{2}.nox];

figure
t = plot(prob);
