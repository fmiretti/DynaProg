%[text] # Use additional outputs
%[text] Additional outputs allow you to save the value for each stage of any additional quantities that you evaluate in the system and cost function. 
%[text] When the optimization problem is solved, any additional output you specified is stored in the problem structure in a cell array. Each cell contains the value of an additional output at each time stage as a numeric vector.
%[text] #### Include an additional output in the system and cost function
%[text] Create the system and cost function with any of the valid signatures, and add any additional output starting from the fourth output.
[x_next, stage_cost, unfeas, addOutput1, .., addOutputN] = sysfun(x, u, w)
%[text] You can specify any number of additional outputs, but they always must be specified starting from the fourth output as the first three outputs are reserved to the updated state variables, the stage cost and unfeasibilities.
%[text] #### Access additional outputs
%[text] When you define additional outputs in the system and cost function, after you run the optimization, the problem structure also contains `AddOutputsProfile,` which is a cell array containing the profiles of the additional outputs.
%[text] 
%[text] 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
