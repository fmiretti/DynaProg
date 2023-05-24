function t = plot(obj)
%DynaProg.plot Plot optimization results
%   plot(prob) plots state variables, control variables and
%   cumulative cost profiles.

% Check the MATLAB version
if verLessThan('matlab','9.7') % aka 2019b
    compatibility_mode = true;
else
    compatibility_mode = false;
end

if isempty(obj.StateProfile)
    error('DynaProg:notSolved', 'This problem structure does not contain any solution. Run the optimization to produce results.')
end
% Set up layout
ncols = lcm(length(obj.N_SV), length(obj.N_CV));
sv_span = ncols/length(obj.N_SV);
cv_span = ncols/length(obj.N_CV);

if compatibility_mode
    subplot(3, ncols, 1);
else
    t = tiledlayout(3, ncols);
end

ax = [];
% Plot SV profiles
for n = 1:length(obj.N_SV)
    if compatibility_mode
        ax(end+1) = subplot(3, ncols, [1 + sv_span*(n-1) sv_span*n]);
    else
        ax(end+1) = nexttile([1 sv_span]); %#ok<*AGROW>
    end
    if isempty(obj.Time)
        plot(obj.StateProfile{n}, 'LineWidth', 1.5)
    else
        plot(obj.Time, obj.StateProfile{n}, 'LineWidth', 1.5)
    end
    title(obj.StateName(n))
    axis tight
end
% Plot CV profiles
for n = 1:length(obj.N_CV)
    if compatibility_mode
        ax(end+1) = subplot(3, ncols, [ncols + 1 + cv_span*(n-1), ncols + cv_span*n]);
    else
        ax(end+1) = nexttile([1 cv_span]); %#ok<*AGROW>
    end
    if isempty(obj.Time)
        plot(obj.ControlProfile{n}, 'LineWidth', 1.5)
    else
        plot(obj.Time(1:end-1), obj.ControlProfile{n}, 'LineWidth', 1.5)
    end
    title(obj.ControlName(n))
    axis tight
end
% Plot cumulative cost profile
if compatibility_mode
    ax(end+1) = subplot(3, ncols, [ncols*2 + 1, ncols*3]);
else
    ax(end+1) = nexttile([1 ncols]); %#ok<*AGROW>
end
if isempty(obj.Time)
    cumCost = cumsum(obj.CostProfile(1:end-1));
    plot(cumCost, 'LineWidth', 1.5)
    xlabel('Stage number')
else
    cumCost = cumsum(obj.CostProfile(1:end-1));
    plot(obj.Time(1:end-1), cumCost, 'LineWidth', 1.5)
    xlabel('Time')
end
title(obj.CostName)
axis tight

% Finalize
if compatibility_mode
    arrayfun(@(x) set(x, 'FontSize', 10), ax)
    t = gcf;
else
    arrayfun(@(x) set(x, 'FontSize', 10), t.Children)
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
end
linkaxes(ax,'x')
end
