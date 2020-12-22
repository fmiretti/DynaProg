function [x_new, stageCost, unfeas] = sysname_int(x, u, w, intVar)
% Basic internal system and cost function template for split models
% Copy this file to your working folder and make sure to rename it.

% State update
dt = 1;
x_new{1} = x{1} .* intVar{1};
x_new{2} = x{2} .* u{2} ./ w{1} + intVar{2};

% Stage cost
stageCost = intVar{3}.*dt;

% unfeasibility
unfeas = [];
end