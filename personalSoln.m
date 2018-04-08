%% Personal HW1
% Jean-Christophe Perrin
clear all;
clc;

%% Given Values
% These values vary for SLS versus cruise conditions. They can be modified
% for different scenarios.

TstaticOut = 218.8; % [k]
PstaticOut = 0.239e5; %[Pa]
Ma = 0.78;
PratioOverall = 32; % pressure ratio of fan and compressors
PratioFan = 1.55; % pressure ratio of the fan
Tturbine = 1450; % [K]
MdotAir = 110; % [kg/s]
BPR = 10;

PrecoveryFactor = 0.998;
PLossRatio = 0.95; % Pressure loss in combustor

%% Outside of the Jet
% The static pressures and temperatures of the air in the reference frame
% of the jet are modified by the speed at which the jet is moving through
% the air. The pressure ratio and temperature ratio are both found in a
% table for any given Mach number and need to be looked up by the user.

disp(Ma);
PratioOut = 1.4947; % from table
TratioOut = 1.1217; % from table

P0Ambient = PstaticOut * PratioOut;
disp(P0Ambient);

%% Inlet to Fan
% There is a slight increase in pressure as the air is diffused before the
% fan.

P0Fan = P0Ambient * PrecoveryFactor;
disp(P0Fan);

%% Input to combustor

P0Combustion = P0Fan * PratioOverall;
disp(P0Combustion);

%% Input to Turbines

P0Turbines = P0Combustion * PLossRatio;
disp(P0Turbines);

%% Input to Nozzle


