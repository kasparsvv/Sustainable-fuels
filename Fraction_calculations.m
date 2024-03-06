clear all;
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
%% Load E5 Mixture
cFuelEthanol = 'C2H5OH';
iSpE5 = myfind({Sp.Name},{cFuelGasoline,cFuelEthanol,'O2','CO2','H2O','N2'});
SpSE5=Sp(iSpE5);                                                              % Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSpE5 = length(SpSE5);
MiE5 = [SpSE5.Mass];
%% Calculation Rg pure gasoline
% C7.76H13.1 + 11.035 O2 + 41.5 N2 --> 7.76 CO2 + 6.55 H2O + 41.5 N2
N_fuel = 1;
N_Oxygen = 11.035;
N_Nitrogen_before_combustion = (11.035*0.79)/0.21;
N_total_before_comb = N_fuel + N_Oxygen + N_Nitrogen_before_combustion;
m_total_before_comb = N_fuel*MiGasoline(1) + N_Oxygen*MiGasoline(2) + N_Nitrogen_before_combustion*MiGasoline(5);
M_total = m_total_before_comb/N_total_before_comb;
Rg_before_comb = Runiv/M_total; %Rg = 275

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
Rg_after_comb = Runiv/M_total_after_comb;
N_frac_CO2 = N_CO2/N_total_after_comb;
N_frac_H20 = N_H20/N_total_after_comb;
N_frac_N2_after_comb = N_Nitrogen_after_comb/N_total_after_comb;
X_comb_out = [0 0 N_frac_CO2 N_frac_H20 N_frac_N2_after_comb];
Y_comb_out = (X_comb_out.*MiGasoline)/(M_total_after_comb);

%% Calculation Rg E5
% C7.76H13.1 + 11.035 O2 + 41.5 N2 --> 7.76 CO2 + 6.55 H2O + 41.5 N2
%C2H5OH + 3 O2 → 2 CO2 + 3 H2O
N_Gasoline_E5 = MolesGasolineE5;
N_Oxygen_E5 = MolesOxygenE5Gasoline + MolesOxygenE5Ethanol;
N_Nitrogen_E5_before_comb = MolesNitrogenE5Gasoline + MolesNitrogenE5Ethanol;
N_Ethanol_E5 = MolesEthanolE5;
N_E5_total_before_comb = N_Gasoline_E5 + N_Oxygen_E5 + N_Nitrogen_E5_before_comb + N_Ethanol_E5;
m_E5_total_before_comb = N_Gasoline_E5*MiE5(1) + N_Oxygen_E5*MiE5(3) + N_Nitrogen_E5_before_comb*MiE5(6) + N_Ethanol_E5*MiE5(2);
M_E5_total = m_E5_total_before_comb/N_E5_total_before_comb;
Rg_E5_before_comb = Runiv/M_E5_total;

%% Calculation Yi for E5 before combustion
N_frac_E5Gasoline = N_Gasoline_E5/N_E5_total_before_comb;
N_frac_E5Oxygen = N_Oxygen_E5/N_E5_total_before_comb;
N_frac_E5Nitrogen = N_Nitrogen_E5_before_comb/N_E5_total_before_comb;
N_frac_E5Ethanol = N_Ethanol_E5/N_E5_total_before_comb;
X_E5_comp_in = [N_frac_E5Gasoline N_frac_E5Ethanol N_frac_E5Oxygen 0 0 N_frac_E5Nitrogen];
Y_E5_comp_in = (X_E5_comp_in.*MiE5)/(M_E5_total);

%% Calculation Yi for E5 after combustion
N_CO2_E5 = 7.76*N_Gasoline_E5 + 2*N_Ethanol_E5;
N_H20_E5 = 6.55*N_Gasoline_E5 + 3*N_Ethanol_E5;
N_Nitrogen_E5_after_comb = N_Nitrogen_E5_before_comb;
N_total_E5_after_comb = N_CO2_E5 + N_H20_E5 + N_Nitrogen_E5_after_comb;
m_total_E5_after_comb = N_CO2_E5*MiE5(4) + N_H20_E5*MiE5(5) + N_Nitrogen_E5_after_comb*MiE5(6);
M_total_E5_after_comb = m_total_E5_after_comb/N_total_E5_after_comb;
Rg_E5_after_comb = Runiv/M_total_E5_after_comb;
N_frac_E5_CO2 = N_CO2_E5/N_total_E5_after_comb;
N_frac_E5_H20 = N_H20_E5/N_total_E5_after_comb;
N_frac_E5_N2_after_comb = N_Nitrogen_E5_after_comb/N_total_E5_after_comb;
X_E5_comb_out = [0 0 0 N_frac_E5_CO2 N_frac_E5_H20 N_frac_E5_N2_after_comb];
Y_E5_comb_out = (X_E5_comb_out.*MiE5)/(M_total_E5_after_comb);

%% Calculation Rg E10
% C7.76H13.1 + 11.035 O2 + 41.5 N2 --> 7.76 CO2 + 6.55 H2O + 41.5 N2
% 2 C8H18 + 25 O2 → 16 CO2 + 18 H2O
N_Gasoline_E10 = MolesGasolineE10;
N_Oxygen_E10 = MolesOxygenE10Gasoline + MolesOxygenE10Ethanol;
N_Nitrogen_E10_before_comb = MolesNitrogenE10Gasoline + MolesNitrogenE10Ethanol;
N_Ethanol_E10 = MolesEthanolE10;
N_E10_total_before_comb = N_Gasoline_E10 + N_Oxygen_E10 + N_Nitrogen_E10_before_comb + N_Ethanol_E10;
m_E10_total_before_comb = N_Gasoline_E10*MiE5(1) + N_Oxygen_E10*MiE5(3) + N_Nitrogen_E10_before_comb*MiE5(6) + N_Ethanol_E10*MiE5(2);
M_E10_total = m_E10_total_before_comb/N_E10_total_before_comb;
Rg_E10_before_comb = Runiv/M_E10_total;

%% Calculation Yi for E10 before combustion
N_frac_E10Gasoline = N_Gasoline_E10/N_E10_total_before_comb;
N_frac_E10Oxygen = N_Oxygen_E10/N_E10_total_before_comb;
N_frac_E10Nitrogen = N_Nitrogen_E10_before_comb/N_E10_total_before_comb;
N_frac_E10Ethanol = N_Ethanol_E10/N_E10_total_before_comb;
X_E10_comp_in = [N_frac_E10Gasoline N_frac_E10Ethanol N_frac_E10Oxygen 0 0 N_frac_E10Nitrogen];
Y_E10_comp_in = (X_E10_comp_in.*MiE5)/(M_E10_total);

%% Calculation Yi for E10 after combustion
N_CO2_E10 = 7.76*N_Gasoline_E10 + 2*N_Ethanol_E10;
N_H20_E10 = 6.55*N_Gasoline_E10 + 3*N_Ethanol_E10;
N_Nitrogen_E10_after_comb = N_Nitrogen_E10_before_comb;
N_total_E10_after_comb = N_CO2_E10 + N_H20_E10 + N_Nitrogen_E10_after_comb;
m_total_E10_after_comb = N_CO2_E10*MiE5(4) + N_H20_E10*MiE5(5) + N_Nitrogen_E10_after_comb*MiE5(6);
M_total_E10_after_comb = m_total_E10_after_comb/N_total_E10_after_comb;
Rg_E10_after_comb = Runiv/M_total_E10_after_comb;
N_frac_E10_CO2 = N_CO2_E10/N_total_E10_after_comb;
N_frac_E10_H20 = N_H20_E10/N_total_E10_after_comb;
N_frac_E10_N2_after_comb = N_Nitrogen_E10_after_comb/N_total_E10_after_comb;
X_E10_comb_out = [0 0 0 N_frac_E10_CO2 N_frac_E10_H20 N_frac_E10_N2_after_comb];
Y_E10_comb_out = (X_E10_comb_out.*MiE5)/(M_total_E10_after_comb);

%% Calculation Rg E15
% C7.76H13.1 + 11.035 O2 + 41.5 N2 --> 7.76 CO2 + 6.55 H2O + 41.5 N2
% 2 C8H18 + 25 O2 → 16 CO2 + 18 H2O
N_Gasoline_E15 = MolesGasolineE15;
N_Oxygen_E15 = MolesOxygenE15Gasoline + MolesOxygenE15Ethanol;
N_Nitrogen_E15_before_comb = MolesNitrogenE15Gasoline + MolesNitrogenE15Ethanol;
N_Ethanol_E15 = MolesEthanolE15;
N_E15_total_before_comb = N_Gasoline_E15 + N_Oxygen_E15 + N_Nitrogen_E15_before_comb + N_Ethanol_E15;
m_E15_total_before_comb = N_Gasoline_E15*MiE5(1) + N_Oxygen_E15*MiE5(3) + N_Nitrogen_E15_before_comb*MiE5(6) + N_Ethanol_E15*MiE5(2);
M_E15_total = m_E15_total_before_comb/N_E15_total_before_comb;
Rg_E15_before_comb = Runiv/M_E15_total;

%% Calculation Yi for E15 before combustion
N_frac_E15Gasoline = N_Gasoline_E15/N_E15_total_before_comb;
N_frac_E15Oxygen = N_Oxygen_E15/N_E15_total_before_comb;
N_frac_E15Nitrogen = N_Nitrogen_E15_before_comb/N_E15_total_before_comb;
N_frac_E15Ethanol = N_Ethanol_E15/N_E15_total_before_comb;
X_E15_comp_in = [N_frac_E15Gasoline N_frac_E15Ethanol N_frac_E15Oxygen 0 0 N_frac_E15Nitrogen];
Y_E15_comp_in = (X_E15_comp_in.*MiE5)/(M_E15_total);

%% Calculation Yi for E15 after combustion
N_CO2_E15 = 7.76*N_Gasoline_E15 + 2*N_Ethanol_E15;
N_H20_E15 = 6.55*N_Gasoline_E15 + 3*N_Ethanol_E15;
N_Nitrogen_E15_after_comb = N_Nitrogen_E15_before_comb;
N_total_E15_after_comb = N_CO2_E15 + N_H20_E15 + N_Nitrogen_E15_after_comb;
m_total_E15_after_comb = N_CO2_E15*MiE5(4) + N_H20_E15*MiE5(5) + N_Nitrogen_E15_after_comb*MiE5(6);
M_total_E15_after_comb = m_total_E15_after_comb/N_total_E15_after_comb;
Rg_E15_after_comb = Runiv/M_total_E15_after_comb;
N_frac_E15_CO2 = N_CO2_E15/N_total_E15_after_comb;
N_frac_E15_H20 = N_H20_E15/N_total_E15_after_comb;
N_frac_E15_N2_after_comb = N_Nitrogen_E15_after_comb/N_total_E15_after_comb;
X_E15_comb_out = [0 0 0 N_frac_E15_CO2 N_frac_E15_H20 N_frac_E15_N2_after_comb];
Y_E15_comb_out = (X_E15_comb_out.*MiE5)/(M_total_E15_after_comb);