function [x_new, stageCost, unfeas] = two_tanks(x, u, ~)
%Two-Tank system example
% Taken from: P. Elbert, S. Ebbesen, and L. Guzzella, “Implementation of
% Dynamic Programming for n-Dimensional Optimal Control Problems With Final
% State Constraints,” IEEE Transactions on Control Systems Technology, vol.
% 21, no. 3, pp. 924–931, May 2013, doi: 10.1109/TCST.2012.2190935.

dt = 0.01; % s

x_new{1} = x{1} + (-0.5*x{1} + u{1}.*u{2})*dt;
x_new{2} = x{2} + (-0.5*x{2} + u{1}.*(1-u{2}))*dt;

% Stage cost
stageCost = (u{1} + 0.1*abs(u{2}-0.5)).*dt;

% unfeasibility
unfeas = [];
end