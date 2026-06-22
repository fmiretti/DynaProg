function [x_next, stageCost] = wrong_outputs_sys(x, u, ~)
% Helper for tInputValidation: declares only 2 outputs (3 are required).

x_next = {x{1}};

stageCost = u{1}.^2;

end
