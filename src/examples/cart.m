function [x_new, stageCost, unfeas] = cart(x, u, ~)
dt = 0.1;
mass = 1;

x_new{1} = x{1} + x{2}.*dt + zeros(size(u{1}));
x_new{2} = x{2} + (1/mass).*u{1}.*dt;

% Stage cost
stageCost = (u{1}.^2).*dt;

% unfeasibility
unfeas = [];
end