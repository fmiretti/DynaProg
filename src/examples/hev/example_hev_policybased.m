%HEV example - policy-based vs. value-based forward algorithm
% Compares the default value-based algorithm with the policy-based one on
% the same p2 parallel HEV problem. 
% The algorithm differs in the forward phase:
%   - Value-based: selects the optimal control from the discrete control 
%       grids.
%   - Policy-based: looks up a pre-computed control map. Control variables
%       marked with the property ControlType = 'continuous' are linearly
%       interpolated and can take any value in range, not just the grid 
%       points.

%% Set up the problem (Value-based) 
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
u2_grid = -1:0.25:1;
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

% Define handle for the sys&cost fuction
sysfun = @(x, u, w) hev(x, u, w, veh, fd, gb, eng, em, batt);

%% Solve with value-based algorithm
prob_vb = DynaProg(x_grid, x_init, x_final, u_grid, Nint, sysfun, ...
    'ExogenousInput', w);

prob_vb = run(prob_vb);

%% Solve with policy-based algorithm
prob_pb = DynaProg(x_grid, x_init, x_final, u_grid, Nint, sysfun, ...
    'ExogenousInput', w, ...
    'ForwardMode',   'policyBased', ...
    'ControlType',   {'discrete', 'continuous'});

prob_pb = run(prob_pb);

%% Compare results
dist_km = trapz(time_s, speed_kmh ./ 3600);
fc2fe = @(grams) (grams / eng.fuel_den) / dist_km * 100; 

FE_vb = prob_vb.totalCost ./ eng.fuel_den / dist_km * 100;
FE_pb = prob_pb.totalCost ./ eng.fuel_den / dist_km * 100;

fprintf('Fuel economy:\n')
fprintf('\tValue-based:  %.2f l/100km\n', FE_vb)
fprintf('\tPolicy-based: %.2f l/100km\n', FE_pb)

%% Plot comparison
figure
tiledlayout(3, 1)

nexttile
plot(time_s, prob_vb.StateProfile{1}(1:end-1), 'k',  'DisplayName', 'Value-based')
hold on
plot(time_s, prob_pb.StateProfile{1}(1:end-1), 'r', 'DisplayName', 'Policy-based')
ylabel('SOC')
legend

nexttile
plot(time_s, prob_vb.ControlProfile{2}, 'k', 'DisplayName', 'Value-based')
hold on
plot(time_s, prob_pb.ControlProfile{2}, 'r', 'DisplayName', 'Policy-based')
ylabel('Torque split')
legend

nexttile
plot(time_s, prob_vb.ControlProfile{1}, 'k', 'DisplayName', 'Value-based')
hold on
plot(time_s, prob_pb.ControlProfile{1}, 'r', 'DisplayName', 'Policy-based')
ylabel('Gear number')
legend
