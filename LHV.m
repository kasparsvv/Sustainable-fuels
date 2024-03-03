%% Basic Code for Retrieving Nasa Tables
relativepath_to_generalfolder='General'; % relative reference to General folder (assumes the folder is in you working folder)
addpath(relativepath_to_generalfolder); 
% Load Nasadatabase
TdataBase=fullfile('General','Nasa','NasaThermalDatabase');
load(TdataBase);

global Runiv Pref
Runiv=8.314472;
Pref=1.01235e5; % Reference pressure, 1 atm
Tref=298.15;    % Reference Temperature
%% Input Values for Calculating Q_LHV
cfuel = 'Gasoline'; %Use Chemical Formula except for 'Gasoline','Diesel'and 'LPG'
T = 20+273.15; %Input Temperature for Q_LHV

%% Extract Specific Nasa Tables, has to be outside the Function
%This is not necessary since Toon already calculated it
%iSp = myfind({Sp.Name},{'O2','CO2','H2O','N2'});     
%Therefore:
run("Fraction_calculations.m");
run("AFcalculations.m");
FuelIndex = myfind({Sp.Name},{cFuelGasoline}); 
%iSp = iSpGasoline
%% Function for Calculating Q_LHV
function [Q_LHV] = LowerHeatingValue(T,SpSGasoline,iSpGasoline)
MoleH2O = 6.55;
MoleCO2 = 7.76;
MoleO2 = 11.035;
MoleN2 = 41.5;
MassH2O = (MoleH2O*0.0180)/SpSGasoline.Mass;
MassCO2 = (MoleCO2*0.0440)/SpSGasoline.Mass;
MassO2 = (MoleO2*0.0320)/SpSGasoline.Mass;
MassN2 = (MoleN2*0.0280)/SpSGasoline.Mass;
Q_LHV = -MassN2*HNasa(T,SpSGasoline(4))-MassCO2*HNasa(T,SpSGasoline(2)) - MassH2O*HNasa(T,SpSGasoline(3)) + MassO2*HNasa(T,SpSGasoline(1)) + HNasa(T,SpSGasoline);
end

