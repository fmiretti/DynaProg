function [x_next, stageCost, unfeas] = trivial_sys(x, u, ~)
% Trivial single model: x_next = u, stage cost = u^2.
% Third argument (exogenous input) is suppressed.
x_next    = {u{1}};
stageCost = u{1}.^2;
unfeas    = [];
