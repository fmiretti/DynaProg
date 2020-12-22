function [intVar, unfeas] = sysname_ext(u, w)
% Basic external system and cost function template for split models
% Copy this file to your working folder and make sure to rename it.

% Evaluate intermediate variables
dt = 1;
a = u{1} + u{2};
b = (u{1} .* u{2}) ./ w{1};

% unfeasibility
unfeas = [];

% Store intermediate variables
intVar{1} = a;
intVar{2} = b;

end