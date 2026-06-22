function [x_next, stageCost, unfeas] = trivial_additive_sys(x, u, ~)
% Trivial additive model: x_next = x + u, stage cost = u^2.
% Third argument (exogenous input) is suppressed.
x_next    = {x{1} + u{1}};
stageCost = u{1}.^2;
unfeas    = [];
