function [m, unfeas] = trivial_ext(u, ~)
% Trivial external model (split-model): m = u (passthrough).
% Second argument (exogenous input) is suppressed.
m      = {u{1}};
unfeas = [];
