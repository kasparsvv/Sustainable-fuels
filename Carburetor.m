%% Load the NASA tables

relativepath_to_generalfolder='General'; % relative reference to General folder (assumes the folder is in you working folder)
addpath(relativepath_to_generalfolder); 
% Load Nasadatabase
TdataBase=fullfile('General','Nasa','NasaThermalDatabase');
load(TdataBase);
run("AFcalculations.m");
run("Fraction_calculations.m");
run("system_parameters.m")

%% Calculate the mass flow

Cd = 0.75;  % OrificeFlowCoefficient, usually between 0.6 and 0.9. Needs to be looked into
Ao = 1;     % [m^2] OrificeSurfaceArea, temporary value
p1 = 1;     %[bar] pressureFloatchamber, 1 atm (not sure if this is the atmospheric pressure)
p2 = 0.5;   %[bar] pressureAirtube, random temporary value, will need a for loop since it changes over time.
densityfuel = MGasoline/1;
 
%mass flow of the mixture using Bernoulli
massflowFuel = Cd * Ao * sqrt(densityfuel * (p1 - p2)); %[kg/s]

%% Calculate AF ratio and Lambda
massflowAir = 1; %[kg/s] actual value changes so needs to be modeled

AF = massflowAir/massflowFuel;
AFstoich = AirFuelRatioGasoline; %value should be in AFcalculations
lambda = AF / AFstoich; % This value for pure gasoline will be close to 1
