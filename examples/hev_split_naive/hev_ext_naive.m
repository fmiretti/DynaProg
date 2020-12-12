function [intVar, unfeas] = hev_ext_naive(u, w, veh, fd, gb, eng, em)
%% Vehicle Model
% Wheels
% Wheel speed (rad/s)
wheelSpd  = w{1} ./ veh.wh_radius;
% Wheel acceleration (rad/s^2)
wheelAcc = w{2} ./ veh.wh_radius;
% Tractive Force (N)
rolling_friction = veh.mass .* veh.gravity .* (veh.first_rrc + veh.second_rrc.*w{1});
vehForce = (wheelSpd~=0) .* (rolling_friction + veh.aero_coeff.*w{1}.^2 + veh.mass.*w{2});
% Wheel torque (Nm)
wheelTrq = (vehForce .* veh.wh_radius + veh.axle_loss .* (wheelSpd~=0));
% Final Drive
% Final drive input speed (rad/s)
fdSpd = fd.spdRatio .* wheelSpd;
% Final drive input acceleration (rad/s^2)
fdAcc = fd.spdRatio .* wheelAcc;
% Final drive input torque (Nm)
fdTrq  = wheelTrq ./ fd.spdRatio ...
    + fd.loss.*(wheelTrq>0) - fd.loss.*(wheelTrq<=0);

% Gearbox
gbSpRatio = gb.spdRatio(u{1});

% Crankshaft speed (rad/s)
shaftSpd  = gbSpRatio .* fdSpd;
% Crankshaft acceleration (rad/s^2)
shaftAcc = gbSpRatio .* fdAcc;
% Gearbox efficiency (-)
gbEff = gb.effMap(fdSpd, fdTrq, u{1});
gbEff = min(max(gbEff, eps), 1);
% Gearbox inertia
gbInertia = gb.inertia .* shaftAcc;
% Crankshaft torque (Nm)
gbLoss = (fdTrq>0 & gbEff>eps)  .* (1-gbEff) .* fdTrq ./ (gbEff .* gbSpRatio) ...
     + (fdTrq>0 & gbEff==eps)  .* fdTrq .* gbSpRatio ...
    + (fdTrq<=0) .* (1-gbEff) .* fdTrq ./ gbSpRatio;

shaftTrq  = fdTrq ./ gbSpRatio + gbLoss + gbInertia;

% Torque Split
% Engine inertia torque (Nm)
engResTrq  = shaftAcc * eng.inertia;
% Electric motor drag torque (Nm)
emResTrq  = shaftAcc * em.inertia;
% Total required torque (Nm)
reqTrq = engResTrq.*(u{2}~=1) + emResTrq + shaftTrq;
% Torque provided by engine
engTrq  = (shaftSpd>0) .* (reqTrq>0)  .* (1-u{2}).*reqTrq;
% Torque provided by electric motor
emTrq  = (shaftSpd>0) .* u{2} .* reqTrq;

pwtUnfeas = (reqTrq<0 & emTrq>0);

% Torque-coupling device (ideal)
emSpd = shaftSpd .* 1.74;
emTrq = emTrq ./ 1.74;

% Engine
% Fuel and pollutants mass flow rates
engSpd = shaftSpd;
fuelFlwRate = eng.fuelMap(engSpd, engTrq); % fuel flow rate, g/s
hcFlwRate = eng.hcMap(engSpd, engTrq); % g/s
coFlwRate = eng.coMap(engSpd, engTrq); % g/s
noxFlwRate = eng.noxMap(engSpd, engTrq); % g/s
egTemp = eng.exhTempMap(engSpd, engTrq); % exhuast gases temperature, C
fuelFlwRate(engTrq==0) = 0;
hcFlwRate(engTrq==0) = 0;
coFlwRate(engTrq==0) = 0;
noxFlwRate(engTrq==0) = 0;
% Maximum engine torque
engMaxTrq = eng.maxTrq(engSpd);
% Costraints
engUnfeas = (engTrq > engMaxTrq) | (engTrq>0 & engSpd < min(eng.spdBrk)) | (engTrq>0 & engSpd > max(eng.spdBrk));

% EM
% Electric motor efficiency
emEff = (shaftSpd~=0) .* em.effMap(emSpd, emTrq) + (shaftSpd==0);
emEff(isnan(emEff)) = 1;
% Calculate electric power consumption
emElPwr =  (emTrq<0) .* emSpd.*emTrq.*emEff + (emTrq>=0) .* emSpd.*emTrq./emEff;
% Limit Torque
emMaxTrq = em.maxTrq(emSpd);
emMinTrq = em.minTrq(emSpd);
% Constraints 
emUnfeas = (isnan(emEff)) + (emTrq<0)  .* (emTrq < emMinTrq) +...
                   (emTrq>=0) .* (emTrq > emMaxTrq);

%% Store intermediate variables
intVar{1} = emElPwr;
intVar{2} = fuelFlwRate;
intVar{3} = egTemp;
intVar{4} = hcFlwRate;
intVar{5} = coFlwRate;
intVar{6} = noxFlwRate;

%% Combine unfeasibilites
unfeas = (pwtUnfeas | engUnfeas | emUnfeas);
