%[text] # Troubleshooting the HEV with a Three-Way Catalyst example
%[text] Consider the [HEV with a Three-Way Catalyst](file:./example_split_hev.mlx). The system and cost functions that are used there are designed with some extra-care so as to avoid all the common issues that may arise with DynaProg. These issues are tipically related to the way in which MATLAB handles n-dimensional arrays; for more information, see [Troubleshooting your model and Safe Mode](file:./troubleshooting.mlx).
%[text] We will start from a naively designed model function, analyze some issues that prevent their usage without using Safe Mode and discuss some solutions.
%%
%[text] ## Troubleshooting the system and cost functions
%[text] Open the naive version of the external model by entering
open('hev_ext_naive.m')
%[text] in the Command Window.
%[text] Since this optimization problem involves two stae variables and two control variable, which are discretized on computational grids having $N\_{x\_1}${"altText":"N\_(x\_1)"}, $N\_{x\_2}${"altText":"N\_(x\_2)"}, $N\_{u\_1}${"altText":"N\_(u\_1)"} and $N\_{u\_2}${"altText":"N\_(u\_2)"} elements respectively, DynaProg will use the external model function passing as inputs:
%[text] - `u{1}` as an array of size $1 \\times 1 \\times N\_{u\_1} \\times 1 ${"altText":"1 x 1 x N\_(u\_1) x 1"},
%[text] - `u{1}` as an array of size $1 \\times 1 \\times 1 \\times N\_{u\_2}${"altText":"1 x 1 x 1 x N\_(u\_2)"}. \
%[text] And it will expect the outputs, `intVar` and `unFeas,` to be arrays of size $1 \\times 1 \\times N\_{u\_1} \\times N\_{u\_2}${"altText":"1 x 1 x N\_(u\_1) x 1"},  $1 \\times 1 \\times N\_{u\_1} \\times 1 ${"altText":"1 x 1 x N\_(u\_1) x 1"} or $1 \\times 1 \\times 1 \\times N\_{u\_2}${"altText":"1 x 1 x 1 x N\_(u\_2)"}.
%[text] ### Array indexing
%[text] The first issue lies in line 23 of the external function.
% Gearbox
gbSpRatio = gb.spdRatio(u{1});
%[text] The issue here is that DynaProg will index into `gb.spdRatio` with `u{1}` which is a $1 \\times 1 \\times N\_{u\_1} \\times 1 ${"altText":"1 x 1 x N\_(u\_1) x 1"} array, but the result (`gbSpRatio`) will be a $ N\_{u\_1} \\times 1 \\times 1 \\times 1${"altText":"1 x 1 x N\_(u\_1) x 1"}array. This is unintended and DynaProg will interpret this as meaning that the gearbox speed ratio is a function of the first state variable only, rather than being a function of the first control variable only and this in turn will produce an error.
%[text] To fix this, we must make sure that `gbSpRatio` is evaluated as a $1 \\times 1 \\times N\_{u\_1} \\times 1 ${"altText":"1 x 1 x N\_(u\_1) x 1"} array. One way to do this is to explicitly initialize it as having to have the same size as `u{1}`, and then using `(:)` in the left-hand size when assigning values to it to ensure we preserve its shape.
% Gearbox
gbSpRatio = zeros(size(u{1}));
gbSpRatio(:) = gb.spdRatio(u{1});
%[text] You can read more about the way MATLAB handles indexing [in this blog post](https://blogs.mathworks.com/loren/2006/08/09/essence-of-indexing/?s_tid=srchtitle).
%[text] From a physical point of view, this line evaluates the speed ratio of the gear engaged by the vehicle's gearbox as a function of the gear number, that is the control variable $u\_1${"altText":"u\_1"}; `gb.spdRatio` is the vector of speed ratios, with each element being the speed ratio associated to a certain gear number. This vector of speed ratios is created in the `hev_data` function and then passed to the model function as an [additional input](file:./use_additional_inputs.mlx).
%[text] ### Querying interpolants
%[text] At line 30 of the external function, we use a `griddedInterpolant` object to evaluate `gbEff` as a function of two previously computed variables, `fdSpd` and `fdTrq`, and the control variable `u{1}`. 
% Gearbox efficiency (-)
gbEff = gb.effMap(fdSpd, fdTrq, u{1});
%[text] Since, based on previous operations, `fdSpd` and `fdTrq` are not dependant on any state or control variable, they are scalars, while  `u{1}` is a $1 \\times 1 \\times N\_{u\_1} \\times 1 ${"altText":"1 x 1 x N\_(u\_1) x 1"} array. The `griddedInterpolant` object requires that all query variables (inputs) have the same size and it will therefore throw an error.
%[text] To fix this, we expand `fdSpd` and `fdTrq` to have the same size as `u{1}.`
% Gearbox efficiency (-)
gbEff = gb.effMap(fdSpd.*ones(size(u{1})), fdTrq.*ones(size(u{1})), u{1});
%[text] From a physical point of view, this line evaluates the gearbox transmission efficiency as a function of its output shaft speed and torque and the engaged gear. This efficiency is evaluated by interpolating an efficiency map which is created in the `hev_data` function and then passed to the model function as an [additional input](file:./use_additional_inputs.mlx).
%[text] Similar issues arise at lines 82 and 78 of the external function.
%[text] ### Logical indexing
%[text] Open the naive version of the internal model by entering
open('hev_int_naive.m')
%[text] At line 60 of the internal model, we use logical indexing.
RR_HC = (intVar{4} - hcTailpipeFlwRate) ./ egFlwRate .* 29 / 44.1;
RR_HC(egFlwRate==0) = 0;
%[text] At lines 57 and 60 of the internal model, we evaluate the molar fractions of the converted polltants in the exhaust gas flow rate. Since the mass flow rate of exhaust gases `egFlwRate` might be zero, the result of the first line may result in `nan` since we might divide by zero. To avoid this, we add a second line where we set the molar fractions to zero when `egFlwRate` is zero by logical indexing. 
%[text] Based on previous operations, the size of `egFlwRate` is $1 \\times 1 \\times N\_{u\_1} \\times N\_{u\_2} ${"altText":"1 x 1 x N\_(u\_1) x N\_(u\_2)"} since it is a function of the control variables only, while the size of `RR_HC` is $1 \\times N\_{x\_2} \\times N\_{u\_1} \\times N\_{u\_2} ${"altText":"1 x N\_(x\_2) x N\_(u\_1) x N\_(u\_2)"} since it is also a function of the second state variable (the TWC temperature). Therefore, logical indexing will not work as intended here and it will yeld unintended results.
%[text] To fix this, we expand `egFlwRate`  to have the same size as `RR_HC` before using logical indexing.
egFlwRate = egFlwRate + zeros(size(RR_HC));
RR_HC(egFlwRate==0) = 0;
%[text] ## Using Safe Mode
%[text] Troubleshooting the system and costs function allows to solve various typical issues that can interfere with DynaProg's algorithm. This procedure is strongly recommended if you have a complex optimization problems with many state and control variables.
%[text] If, however, your optimization problem is relatively simple or computational time is not critical to your task, you can enable Safe Mode. Safe Mode attempts to automatically handle all these issues at the cost of increased optimization time.
%[text] To see an example, type
open('hev_int_naive.m')
%[text] in the Command Window. 
%[text] 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":15.9}
%---
