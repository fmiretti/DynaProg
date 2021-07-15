%HEV example #1
% Control strategy optimization of a p2 parallel HEV passenger car
% Data source: ADVISOR software (2003, "ADVISOR". National Renewable Energy
% Laboratory. Available: https://sourceforge.net/projects/adv-vehicle-sim/,
% accessed on: Nov 2020)

%% Set up the problem
clear
% State variable grid
SVnames = 'SOC';
x_grid = {0.4:0.001:0.7};
% Initial state 
x_init = {0.6};
% Final state constraints
x_final = {[0.599 0.601]};

% Control variable grid
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
[veh, fd, gb, eng, em, batt] = hev_data();

% Create DynaProg object
prob = DynaProg(x_grid, x_init, x_final, u_grid, Nint, ...
 @(x, u, w) hev(x, u, w, veh, fd, gb, eng, em, batt), 'ExogenousInput', w);

%% Solve and visualize results
% Solve the problem
prob = run(prob);

% Add time vector to use in the plot
prob.Time = [time_s; time_s(end)+dt];
% Add SV, CV and cost names to be used in the plot
prob.StateName = SVnames;
prob.ControlName = CVnames;
prob.CostName = 'Fuel Consumption, g';

% Retrieve our additional outputs
engTrq = prob.AddOutputsProfile{1};
emTrq = prob.AddOutputsProfile{2};

% Plot results
figure
t = plot(prob);

% Edit the plot (only if your MATLAB version supports 'tiledlayout')
if ~verLessThan('matlab','9.6')
    nexttile(t, 5)
    plot(time_s, engTrq, 'LineWidth', 1.5)
    hold on
    plot(time_s, emTrq, 'LineWidth', 1.5)
    legend('Engine torque, Nm', 'EM torque, Nm', 'FontSize', 10)
end
