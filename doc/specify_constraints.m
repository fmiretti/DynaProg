%[text] # Specify contraints on the state and control variables
%[text] In order to specify constraints on the state and control variables, you can introduce unfeasibilty in your system function.
[x_next, stage_cost, unfeas] = sysfun(x, u, ~)
%[text] Unfeasibility, which must be the third output of your function, is a logical variable, and it can be any function of the state variables and control variables.
%[text] In your function, use logical indexing to set it to `true` for those conditions which must be unfeasible, according to your constraints.
%[text] **Examples**
%[text] Enforce that the value of $x\_1${"altText":"x\_1"} does not fall below `0`.
unfeas(x{1}<0) = true;
%[text] Enforce that the value of $u\_1${"altText":"u\_1"} is not `1` if $x\_1${"altText":"x\_1"} is below `0.5`.
unfeas(x{1}<0 & u{1}==1) = true;
%[text] 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
