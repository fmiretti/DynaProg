function [x_new, stageCost, unfeas, engTrq, emTrq] = hev(x, u, w, veh, fd, gb, eng, em, batt)
dt = 1; % s

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
brakeTrq  = (shaftSpd>0) .* (reqTrq<=0) .* (1-u{2}).*reqTrq;
% Torque provided by electric motor
emTrq  = (shaftSpd>0) .* u{2} .* reqTrq;

pwtUnfeas = (reqTrq<0 & emTrq>0);

% Torque-coupling device (ideal)
emSpd = shaftSpd .* 1.74;
emTrq = emTrq ./ 1.74;

% Engine
% Fuel and pollutants mass flow rates
engSpd = shaftSpd.*ones(size(engTrq));
fuelFlwRate = eng.fuelMap(engSpd, engTrq); % fuel flow rate, g/s
fuelFlwRate(engTrq==0) = 0;

% Maximum engine torque
engMaxTrq = eng.maxTrq(engSpd);
% Costraints
engUnfeas = (engTrq > engMaxTrq) | (engTrq>0 & engSpd < min(eng.spdBrk)) | (engTrq>0 & engSpd > max(eng.spdBrk));

% EM
% Electric motor efficiency
emSpd = emSpd.*ones(size(emTrq));
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

% BATTERY
% Assume 700 W of electrical accessories
battPwr = emElPwr + 700;
% Assume inverter efficiency of 0.95
battPwr = (battPwr>0)  .* battPwr ./ 0.95...
  + (battPwr<=0) .* battPwr .* 0.95;

% columbic efficiency
battColumbicEff = (battPwr>0) + (battPwr<=0) .* batt.coulombic_eff;
% Battery internal resistance
battR = (battPwr>0)  .* batt.dischrgRes(x{1})...
  + (battPwr<=0) .* batt.chrgRes(x{1});
% Battery voltage
battVoltage = batt.ocv(x{1});

% Battery current
battCurr = battColumbicEff .* (battVoltage-sqrt(battVoltage.^2 - 4.*battR.*battPwr))./(2.*battR);
battCurr = real(battCurr);
% Maximum charge current
maxChrgBattCurr = (battPwr<=0) .* batt.minCurr(x{1});
% Maximum discharge power
maxBattDisPwr = (battPwr>0) .* batt.maxPwr(x{1});

% New battery state of charge
x_new{1}  = - battCurr ./ (batt.maxCap * 3600) .* dt + x{1};
% Constraints
battUnfeas = (battPwr<=0).*(battCurr < maxChrgBattCurr) + (battPwr>0).*(battPwr > maxBattDisPwr);

%% Stage cost
stageCost  = fuelFlwRate;

%% Combine unfeasibilites
unfeas = (pwtUnfeas | engUnfeas | emUnfeas | battUnfeas);

