function [veh, fd, gb, eng, em, batt, twc] = hev_data()
%% Engine data

% FC_SI41_emis.m - 'Geo 1.0L (41kW) SI Engine - transient data' 
eng.disp = 1.0;  % (L), engine displacement
% (rad/s), speed range of the engine
eng.spdBrk = [104.5 149.2 220.9 292.5 364.1 435.7 507.4 552.2 596.9]; 
% (N*m), torque range of the engine
eng.trqBrk = [6.8 13.6 20.4 27.2 33.8 40.6 47.4 54.2 61 67.8 74.6 81.4]; 

% (g/s), fuel use map indexed vertically by eng.spdBrk and 
% horizontally by eng.trqBrk
eng.fuelMap_gpkWh =[
635.7	635.7	541.4	447.2	352.9	332.2	311.4	322.4	333.5	333.5	333.5	333.5
678.4	500.1	443.8	387.4	331.1	301.8	297	283.4	269.8	358	358	358
463.4	463.4	407.6	350.1	294.3	280.8	267.3	253.9	269.8	303.2	336.7	336.7
699.1	567.9	500.3	432.7	301.4	283.9	266.3	248.7	258.8	268.8	271.9	317.9
592.9	592.9	494.6	393.4	295.1	279.4	263.6	247.9	255.2	262.5	295	322.6
667.9	524.8	381.6	351.9	322.2	304.9	287.5	270.8	290.8	310.9	330.9	330.9
630.6	630.6	522.5	411.1	303	304.4	305.8	304.2	314.5	324.8	327.7	327.7
698.4	500.5	428.6	392.7	356.8	337.9	328.4	319	328.8	338.6	333.7	333.7
751.1	637.8	521.1	407.8	393.1	378.4	363.3	348.2	318.8	340.2	340.2	340.2];
% fuel map in g/kWh 

% (g/s), engine out HC emissions indexed vertically by eng.map_spd and
% horizontally by eng.map_trq
eng.hcMap_gpkWh =[
11.5	11.5	9.8	8.2	6.5	5.8	5.1	5.9	6.8	6.8	6.8	6.8
7.8	7	6.2	5.5	4.7	4.3	4.7	4.6	4.5	4.6	4.6	4.6
5.8	5.8	5.2	4.6	4	4	4	4	4.5	4.6	4.6	4.6
7.1	6	5.4	4.9	3.8	3.7	3.6	3.4	3.2	3	3.4	3.9
5.8	5.8	5	4.3	3.6	3.6	3.6	3.7	3.6	3.6	4	3.9
5.6	4.7	3.7	3.7	3.7	3.4	3.1	3	3.2	3.4	3.5	3.5
8.2	8.2	6.8	5.4	4.1	3.7	3.3	3.1	3.2	3.2	3.3	3.3
5.8	5.2	4.8	4.5	4.3	3.7	3.4	3.2	3.3	3.3	3.3	3.3
5.6	5.8	5.9	6.1	5.7	5.4	5	4.6	3.9	3.9	3.9	3.9];
% engine out HC in g/kWh

% (g/s), engine out CO emissions indexed vertically by eng.map_spd and
% horizontally by eng.map_trq
eng.coMap_gpkWh =[
71.8	71.8	58.8	45.7	32.7	27.5	22.4	82.9	143.3	143.3	143.3	143.3
104.4	68.3	56.8	45.3	33.9	23.3	25.7	24.2	22.8	268.6	268.6	268.6
48.3	48.3	42.9	37.2	31.8	28.6	25.4	22.2	22.8	152.9	283	283
103.1	82.7	72.2	61.8	41.5	36.9	32.3	27.8	31.1	34.4	178.5	279.9
88.1	88.1	74.2	59.9	46	41.9	37.8	33.7	34.8	35.8	158.8	264.6
96.1	74.5	52.9	51.2	49.5	45.9	42.4	34	117.9	201.8	285.7	285.7
114.6	114.6	92.1	69	46.5	60.7	74.8	129.6	195.5	261.4	277.8	277.8
60.1	63.8	64	64.1	64.3	108.2	130.2	152.2	216.7	281.2	278.1	278.1
51.8	75.2	99.3	122.8	134.9	147.1	159.7	172.2	196.6	286.6	286.6	286.6];
% engine out CO in g/kWh

% (g/s), engine out NOx emissions indexed vertically by eng.map_spd and
% horizontally by eng.map_trq
eng.noxMap_gpkWh =[
5.8	5.8	9.3	12.8	16.3	16.1	15.9	13	10.2	10.2	10.2	10.2
5.2	8.8	9.2	9.7	10.2	8	8.9	13.2	17.5	4.6	4.6	4.6
8.1	8.1	8.8	9.6	10.4	10.8	11.3	11.7	17.5	11.6	5.7	5.7
4.2	5.6	6.3	7	8.4	8.9	9.5	10.1	13.9	17.7	8.1	3.1
5.8	5.8	7.2	8.7	10.1	11	11.8	12.6	15.9	19.2	9.3	6.8
14.9	16.4	17.8	19.4	21	20.6	20.3	19.1	14.6	10.2	5.7	5.7
28.7	28.7	26.8	25	23.1	22.4	21.7	16.5	12.1	7.8	6.5	6.5
31.1	27.9	26.7	26.2	25.6	20.9	18.6	16.3	12.1	7.8	6.8	6.8
35	31.1	27.1	23.2	21.1	19.1	17	14.9	10.9	7.4	7.4	7.4];
% engine out NOx in g/kWh

% convert g/kWh maps to g/s maps
[T, spd] = meshgrid(eng.trqBrk, eng.spdBrk);
eng.pwrMap_kW = T.*spd/1000;
eng.fuelMapData = eng.fuelMap_gpkWh.*eng.pwrMap_kW/3600;
eng.hcMapData = eng.hcMap_gpkWh.*eng.pwrMap_kW/3600;
eng.coMapData = eng.coMap_gpkWh.*eng.pwrMap_kW/3600;
eng.noxMapData = eng.noxMap_gpkWh.*eng.pwrMap_kW/3600;
eng.fuelMap = griddedInterpolant({eng.spdBrk, eng.trqBrk}, eng.fuelMapData, 'linear', 'nearest');
eng.hcMap = griddedInterpolant({eng.spdBrk, eng.trqBrk}, eng.hcMapData, 'linear', 'nearest');
eng.coMap = griddedInterpolant({eng.spdBrk, eng.trqBrk}, eng.coMapData, 'linear', 'nearest');
eng.noxMap = griddedInterpolant({eng.spdBrk, eng.trqBrk}, eng.noxMapData, 'linear', 'nearest');

% Torque limit
% (N*m), max torque curve of the engine indexed by eng.map_spd
eng.maxTrqData = [61 67.6 73.7 78.5 80.9 77.3 76.2 73.3 68.7]; 
% (N*m), closed throttle torque of the engine (max torque that can be absorbed)
% indexed by eng.spdBrk -- correlation from JDMA
eng.motTrqData = 4.448/3.281*(-eng.disp)*61.02/24 * ...
   (9*(eng.spdBrk/max(eng.spdBrk)).^2 + 14 * (eng.spdBrk/max(eng.spdBrk)));
eng.maxTrq = griddedInterpolant(eng.spdBrk, eng.maxTrqData);
eng.motTrq = griddedInterpolant(eng.spdBrk, eng.motTrqData);

% other data
eng.inertia = 0.1; % (kg*m^2), rotational inertia of the engine (unknown)
eng.maxPwr = (max(eng.spdBrk.*eng.maxTrqData)/1000); % kW     peak engine power
eng.base_mass = 1.8.*eng.maxPwr;            % (kg), mass of the engine block and head (base engine)
eng.acc_mass = 0.8*eng.maxPwr;             % kg    engine accy's, electrics, cntrl's - assumes mass penalty of 0.8 kg/kW (from OTA report)
eng.fuel_mass = 0.6*eng.maxPwr;            % kg    mass of fuel and fuel tank (from OTA report)
eng.mass = eng.base_mass+eng.acc_mass+eng.fuel_mass; % kg  total engine/fuel system mass

eng.ext_sarea = 0.3*(eng.maxPwr/100)^0.67;       % m^2    exterior surface area of engine
eng.fuel_den = 0.749*1000; % (g/l), density of the fuel 
eng.fuel_lhv = 42.6*1000; % (J/g), lower heating value of the fuel

% calc "predicted" exh gas flow rate and engine-out (EO) temp
eng.ex_pwr_frac = [0.40 0.30];                        % --   frac of waste heat that goes to exhaust as func of engine speed
exflow_mapData = eng.fuelMapData*(1+14.5);                % g/s  ex gas flow map:  for SI engines, exflow = (fuel use)*[1 + (stoic A/F ratio)]
eng.waste_pwr_map = eng.fuelMapData*eng.fuel_lhv - T.*spd;   % W    tot FC waste heat = (fuel pwr) - (mech out pwr)
spd = eng.spdBrk;
eng.ex_pwr_map = zeros(size(eng.waste_pwr_map));       % W   initialize size of ex pwr map
for i = 1:length(spd)
 eng.ex_pwr_map(i,:) = eng.waste_pwr_map(i,:)*interp1([min(spd) max(spd)],eng.ex_pwr_frac,spd(i)); % W  trq-spd map of waste heat to exh 
end
eng.exhTempMapData = eng.ex_pwr_map./(exflow_mapData*1089/1000) + 20 - 273;  % C   EO ex gas temp = Q/(MF*cp) + Tamb (assumes engine tested ~20 C)
eng.exhTempMap = griddedInterpolant({eng.spdBrk, eng.trqBrk}, eng.exhTempMapData);

%% EM data
% MC_AC75 'Westinghouse 75-kW (continuous) AC induction motor/inverter'

% SPEED & TORQUE RANGES over which data is defined
% (rad/s), speed range of the motor
em.spdBrk = [0 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000]*(2*pi/60);
% Conversion from rpm to rad/s 

% (N*m), torque range of the motor
em.trqBrk = [-200 -180 -160 -140 -120 -100 -80 -60 -40 -20 ...
	0 20 40 60 80 100 120 140 160 180 200]*4.448/3.281;

% EFFICIENCY AND INPUT POWER MAPS
% (--), efficiency map indexed vertically by em.map_spd and 
% horizontally by em.trqBrk
em.effMapData = [...
0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7	0.7
0.78	0.78	0.79	0.8	0.81	0.82	0.82	0.82	0.81	0.77	0.7	0.77	0.81	0.82	0.82	0.82	0.81	0.8	0.79	0.78	0.78
0.85	0.86	0.86	0.86	0.87	0.88	0.87	0.86	0.85	0.82	0.7	0.82	0.85	0.86	0.87	0.88	0.87	0.86	0.86	0.86	0.85
0.86	0.87	0.88	0.89	0.9	0.9	0.9	0.9	0.89	0.87	0.7	0.87	0.89	0.9	0.9	0.9	0.9	0.89	0.88	0.87	0.86
0.81	0.82	0.85	0.87	0.88	0.9	0.91	0.91	0.91	0.88	0.7	0.88	0.91	0.91	0.91	0.9	0.88	0.87	0.85	0.82	0.81
0.82	0.82	0.82	0.82	0.85	0.87	0.9	0.91	0.91	0.89	0.7	0.89	0.91	0.91	0.9	0.87	0.85	0.82	0.82	0.82	0.82
0.79	0.79	0.79	0.78	0.79	0.82	0.86	0.9	0.91	0.9	0.7	0.9	0.91	0.9	0.86	0.82	0.79	0.78	0.79	0.79	0.79
0.78	0.78	0.78	0.78	0.78	0.78	0.8	0.88	0.91	0.91	0.7	0.91	0.91	0.88	0.8	0.78	0.78	0.78	0.78	0.78	0.78
0.78	0.78	0.78	0.78	0.78	0.78	0.78	0.8	0.9	0.92	0.7	0.92	0.9	0.8	0.78	0.78	0.78	0.78	0.78	0.78	0.78
0.78	0.78	0.78	0.78	0.78	0.78	0.78	0.78	0.88	0.92	0.7	0.92	0.88	0.78	0.78	0.78	0.78	0.78	0.78	0.78	0.78
0.78	0.78	0.78	0.78	0.78	0.78	0.78	0.78	0.8	0.92	0.7	0.92	0.8	0.78	0.78	0.78	0.78	0.78	0.78	0.78	0.78];
em.effMap = griddedInterpolant({em.spdBrk, em.trqBrk}, em.effMapData);

% LIMITS
% max continuous torque curve of the motor indexed by em.spdBrk
em.maxTrqData = [200 200 200 175.2 131.4 105.1 87.6 75.1 65.7 58.4 52.5]*...
	4.448/3.281; % (N*m)
em.minTrqData = -1*[200 200 200 175.2 131.4 105.1 87.6 75.1 65.7 58.4 52.5]*...
	4.448/3.281; % (N*m)
em.maxTrq = griddedInterpolant(em.spdBrk, em.maxTrqData);
em.minTrq = griddedInterpolant(em.spdBrk, em.minTrqData);

em.inertia = 0;  % (kg*m^2), rotor's rotational inertia
em.mass = 91;  % (kg), mass of motor and controller	

%% Transmission
% TX_5SPD_SI - 'manual 5 speed';
gb.spdRatio = [3.25  1.81  1.21  0.86  0.64]*4.06; % ~'92 Saturn 88-hp 1.9-L SOHC

%TX_VW % FILE ID, LOSSES
load tx_92saturn_26_50_5 % load the tx_eff_map, tx_map_spd, and tx_map_scale into workspace
gb.effMapData = tx_eff_map;
gb.mapSpdBrk = tx_map_spd;
gb.mapTrqBrk = tx_map_trq;
gb.effMap = griddedInterpolant({gb.mapSpdBrk, gb.mapTrqBrk, 1:5}, gb.effMapData);

gb.mass = 141/2.205; % (kg), mass of the gearbox - 1990 Taurus, OTA Report
gb.inertia = 0;	% (kg*m^2), gearbox rotational inertia measured at input; unknown

%final drive variables
fd.loss = 0;    % (Nm), constant torque loss in final drive, measured at input
fd.spdRatio = 1;   % (--), =(final drive input speed)/(f.d. output speed)
fd.inertia = 0; % (kg*m^2), rotational inertia of final drive, measured at input
fd.mass = 110/2.205; % (kg), mass of the final drive - 1990 Taurus, OTA report

%% Battery
% ESS_PB25 - Hawker Genesis 12V26Ah10EP VRLA battery, tested by VA Tech';

% other data
batt.module_mass = 11;  % (kg), mass of a single ~12 V module
batt.module_num = 25;  % a default value for number of modules
batt.mass = batt.module_mass * batt.module_num;

% SOC RANGE over which data is defined
batt.socBrk = [0:.1:1];  % (--)

% LOSS AND EFFICIENCY parameters
% Parameters vary by SOC horizontally, and temperature vertically
batt.maxCap = 25;	% (A*h), max. capacity at C/5 rate, indexed by batt.tempBrk

% average coulombic (a.k.a. amp-hour) efficiency 
batt.coulombic_eff = .9;  % (--);

% module's resistance to being discharged, indexed by batt.soc
batt.dischrgResData = [40.7 37.0 33.8 26.9 19.3 15.1 13.1 12.3 11.7 11.8 12.2]/1000; % (ohm)
batt.dischrgResData = batt.dischrgResData .* batt.module_num;
% module's resistance to being charged, indexed by batt.soc
batt.chrgResData = [31.6 29.8 29.5 28.7 28.0 26.9 23.1 25.0 26.1 28.8 47.2]/1000; % (ohm)
batt.chrgResData = batt.chrgResData .* batt.module_num;
% module's open-circuit (a.k.a. no-load) voltage, indexed by batt.soc
batt.ocvData = [11.70 11.85 11.96 12.11 12.26 12.37 12.48 12.59 12.67 12.78 12.89]; % (V)
batt.ocvData = batt.ocvData .* batt.module_num;
batt.dischrgRes = griddedInterpolant(batt.socBrk, batt.dischrgResData);
batt.chrgRes = griddedInterpolant(batt.socBrk, batt.chrgResData);
batt.ocv = griddedInterpolant(batt.socBrk, batt.ocvData);

% LIMITS
min_volts = 9.5;
max_volts = 16.5;

minCurrData = ( batt.ocvData - max_volts*batt.module_num ) ./ batt.chrgResData;
maxPwrData = batt.ocvData.^2 ./ (4.*batt.dischrgResData);

batt.minCurr = griddedInterpolant(batt.socBrk, minCurrData);
batt.maxPwr = griddedInterpolant(batt.socBrk, maxPwrData);

%% TWC
% EX_SI - 'Standard catalyst for ~stoichiometric SI engine';

% Pollutants conversion efficiencies
twc.tempBrk = [-40 0 220 240 310 415 475 550 650 1200]; % (deg. C)
% catalyst's temperature-dependent conversion efficiencies indexed by ex_cat_temp_range
twc.hcData = [0 0 0 0 0.6 0.84 0.92 1 1 1]*0.85;  % (--)
twc.coData = [0 0 0 0.48 0.82 0.93 0.96 1 1 1]*0.95; % (--)
twc.noxData = [0 0 0.09 0.2 0.89 0.99 0.99 1 1 1]*0.91; % (--)

twc.hcConvEff = griddedInterpolant(twc.tempBrk, twc.hcData);
twc.coConvEff = griddedInterpolant(twc.tempBrk, twc.coData);
twc.noxConvEff = griddedInterpolant(twc.tempBrk, twc.noxData);

% CONVENTIONAL CONVERTER
exhaust_system_mass = 0.26*eng.maxPwr; % kg     mass of exhaust system assumes mass penalty of 0.26 kg/kW (from 1994 OTA report, Table 3)
catalyst_mass = exhaust_system_mass*0.36;           % kg     mass of catalytic converter (from 1994 OTA report, Table 3)
monolith_mass = catalyst_mass*0.22;   % kg     mass of cat monolith (ceramic)
internal_catalyst_mass = catalyst_mass*0.33;   % kg     mass of cat internal SS shell

monolith_cp = 1070;        % J/kgK  ave cp of cat int: CERAMIC MON. (SAE #880282)
internal_catalyst_cp = 460;         % J/kgK  ave cp of cat int: METAL MON. (SAE #890798)

twc.gasSpecHeatCap = 1089;            % J/kgK  ave sens heat cap of exh gas (SAE #890798)

monolith_surface_area = 0.1*(eng.maxPwr/100)^0.67; % m^2  outer surface area of cat monolith (approx. 0.1 m^2/100 kW)
twc.externalArea = monolith_surface_area*1.4;   % m^2  surface area of cat ext shield
twc.externalEmissivity = 0.7;       %        emissivity of cat ext shield
twc.SBconst = 5.670367*1e-8; % Stefan-Boltzmann constant, W m^-2 K^-4
twc.externalConvectionCoeff = 15; % W m^-2 K^-1

twc.heatCapacity = monolith_cp*monolith_mass + internal_catalyst_cp*internal_catalyst_mass;

%% Vehicle Data
% VEH_SMCAR - 'Hypothetical small car' from ADVISOR
veh.gravity=9.81;    % m/s^2
veh.air_density=1.2; % kg/m^3

veh.glider_mass = (2325/2.205)-462; % (kg), vehicle mass w/o propulsion system
veh.cargo_mass = 136; %kg  cargo mass
veh.mass = veh.glider_mass + veh.cargo_mass + eng.mass + em.mass + batt.mass + gb.mass + fd.mass;

veh.CD=0.335;  % (--), coefficient of aerodynamic drag
veh.FA=2.0;    % (m^2), frontal area
veh.aero_coeff = 0.5 .* veh.air_density .* veh.CD .* veh.FA;

% for the eq'n:  rolling_drag=mass*gravity*(veh_1st_rrc+veh_2nd_rrc*v)
veh.first_rrc=0.009;  % (--)
veh.second_rrc=0;		% (s/m)
% wheels
veh.wh_radius=0.282;    % (m), rolling radius of 185/70R14 tire from a ~1990 Mazda 626
veh.wh_inertia=181/2.205*veh.wh_radius^2/2;  % (kg*m^2) 
% drag torque applied at the front (drive) axle
wh_axle_loss_mass = [0 2000];   % (kg)
wh_axle_loss_trq = [4 24]*.4;   % (Nm)
veh.axle_loss = interp1(wh_axle_loss_mass, wh_axle_loss_trq, veh.mass);
