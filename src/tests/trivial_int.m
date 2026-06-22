function [x_next, stageCost, unfeas] = trivial_int(x, u, ~, m)
% Trivial internal model (split-model): x_next = m, stage cost = m^2.
% Third argument (exogenous input) is suppressed.
% Must be called with ExogenousInput so model_wrapper passes 4 args.
x_next    = {m{1}};
stageCost = m{1}.^2;
unfeas    = [];
