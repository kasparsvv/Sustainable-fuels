%% To make sure that matlab will find the functions. You must change it to your situation 
relativepath_to_generalfolder='General'; % relative reference to General folder (assumes the folder is in you working folder)
addpath(relativepath_to_generalfolder); 
run("AFcalculations.m")

%% Load Nasadatabase
TdataBase=fullfile('General','Nasa','NasaThermalDatabase');
load(TdataBase);
global Runiv
Runiv=8.314472;
%% Load Gasoline Mixture
% Composition Gasoline
cFuelGasoline = 'Gasoline';
iSpGasoline = myfind({Sp.Name},{cFuelGasoline,'O2','CO2','H2O','N2'});                    % Find indexes of these species
SpSGasoline=Sp(iSpGasoline);  
NSpGasoline = length(SpSGasoline);
MiGasoline = [SpSGasoline.Mass];
YfuelGasoline = [1 0 0 0 0];
MGasoline = YfuelGasoline*MiGasoline'; 
Gasoline_comp = SpSGasoline.Elcomp; %Answer in order [O H C N AR]
% Result gives: [0 13.1 7.76 0 0]
%% Calculation Rg pure gasoline
% C7.76H13.1 + 11.035 O2 + 41.5 N2 --> 7.76 CO2 + 6.55 H2O + 41.5 N2
N_fuel = 1;
N_Oxygen = 11.035;
N_Nitrogen_before_combustion = (11.035*0.79)/0.21;
N_total_before_comb = N_fuel + N_Oxygen + N_Nitrogen_before_combustion;
m_total_before_comb = N_fuel*MiGasoline(1) + N_Oxygen*MiGasoline(2) + N_Nitrogen_before_combustion*MiGasoline(5);
M_total = m_total_before_comb/N_total_before_comb;
Rg = Runiv/M_total; %Rg = 275

%% Calculation Yi before combustion for Cv
N_frac_fuel = N_fuel/N_total_before_comb;
N_frac_Oxygen = N_Oxygen/N_total_before_comb;
N_frac_Nitrogen = N_Nitrogen_before_combustion/N_total_before_comb;
X_comp_in = [N_frac_fuel N_frac_Oxygen 0 0 N_frac_Nitrogen];
Y_comp_in = (X_comp_in.*MiGasoline)/(M_total);

%% Calculation Yi after combustion for Cp
N_CO2 = 7.76;
N_H20 = 6.55;
N_Nitrogen_after_comb = N_Nitrogen_before_combustion;
N_total_after_comb = N_CO2 + N_H20 + N_Nitrogen_after_comb;
m_total_after_comb = N_CO2*MiGasoline(3) + N_H20*MiGasoline(4) + N_Nitrogen_after_comb*MiGasoline(5);
M_total_after_comb = m_total_after_comb/N_total_after_comb;
Rg = Runiv/M_total_after_comb;
N_frac_CO2 = N_CO2/N_total_after_comb;
N_frac_H20 = N_H20/N_total_after_comb;
N_frac_N2_after_comb = N_Nitrogen_after_comb/N_total_after_comb;
X_comb_out = [0 0 N_frac_CO2 N_frac_H20 N_frac_N2_after_comb];
Y_comb_out = (X_comb_out.*MiGasoline)/(M_total_after_comb);


