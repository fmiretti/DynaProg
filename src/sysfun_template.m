function [x_new, stageCost, unfeas] = sysname(x, u, w)
% Basic system and cost function template
% Copy this file to your working folder and make sure to rename it.

% State update
dt = 1;
x_new{1} = x{1} .* u{1};
x_new{2} = x{2} .* u{2} ./ w{1};

% Stage cost
stageCost = u{1}.*dt;

% unfeasibility
unfeas = [];
end