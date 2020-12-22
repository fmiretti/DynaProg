function [x_new, stageCost, unfeas, fuelFlwRate, tpEmissions] = hev_int_naive(x, ~, ~, intVar, batt, twc)
% x{1}  SOC
% x{2}  TWC temperature C
% intVar{1}  emElPwr W
% intVar{2}  fuelFlwRate g/s
% intVar{3}  egTemp C
% intVar{4}  engine-out hcFlwRate g/s
% intVar{5}  engine-out coFlwRate g/s
% intVar{6}  engine-out noxFlwRate g/s

dt = 1;

% Battery

emElPwr = intVar{1};
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

% Three-way catalyst

% Talipipe pollutants
hcTailpipeFlwRate = (1-twc.hcConvEff(x{2})) .* intVar{4};
coTailpipeFlwRate = (1-twc.coConvEff(x{2})) .* intVar{5};
noxTailpipeFlwRate = (1-twc.noxConvEff(x{2})) .* intVar{6};

fuelFlwRate = intVar{2};
egFlwRate = fuelFlwRate.*(1+14.5); % Exhaust gases flow rate; assumes stoich fuel-air mixture
% Catalist - exhaust gases convective heat transfer
Qconv = egFlwRate .* 1e-3 .* twc.gasSpecHeatCap .* (intVar{3} - x{2});
% Catalist - environment convective heat transfer
Qenv =  - twc.externalConvectionCoeff * twc.externalArea * (x{2} - 20);
% Chem heat release
RR_HC = (intVar{4} - hcTailpipeFlwRate) ./ egFlwRate .* 29 / 44.1;
RR_CO = (intVar{5} - coTailpipeFlwRate) ./ egFlwRate .* 29 / 28;
RR_NOx = (intVar{6} - noxTailpipeFlwRate) ./ egFlwRate .* 29 / 46;
RR_HC(egFlwRate==0) = 0;
RR_CO(egFlwRate==0) = 0;
RR_NOx(egFlwRate==0) = 0;
RR_HC = RR_HC.*100;
RR_CO = RR_CO.*100;
RR_NOx = RR_NOx.*100;
B = 669.*(RR_HC) + (93+79/3).*(RR_CO - RR_NOx) + 121.*RR_NOx;
Qchem = egFlwRate .* 1e-3 .* twc.gasSpecHeatCap .* B;
% Catalist - environment radiation heat transfer
Qrad = - twc.externalArea .* twc.externalEmissivity .* twc.SBconst .* (x{2}.^4 - 20^4);
% Total heat transfer
Qtot = Qconv + Qchem + Qrad + Qenv;
% TWC temperature update
x_new{2} = x{2} + dt .* (1/twc.heatCapacity) .* Qtot;

% Unfeas
unfeas = battUnfeas;

% Cost
stageCost  = fuelFlwRate + (hcTailpipeFlwRate.*10 + coTailpipeFlwRate.*100 + noxTailpipeFlwRate.*100);

% Additional output
tpEmissions.hc = hcTailpipeFlwRate;
tpEmissions.co = coTailpipeFlwRate;
tpEmissions.nox = noxTailpipeFlwRate;
