relativepath_to_generalfolder='General'; % relative reference to General folder (assumes the folder is in you working folder)
addpath(relativepath_to_generalfolder); 
% Load Nasadatabase
TdataBase=fullfile('General','Nasa','NasaThermalDatabase');
load(TdataBase);

DensityGasoline = 705; %kg/m^3
DensityEthanol = 789; %kg/m^3 according to Discussion Page(690-720)
DensityAir = 1.293; %kg/m^3

% Composition Ethanol
cFuelEthanol = 'C2H5OH';        %Ethanol
iSpEthanol = myfind({Sp.Name},{cFuelEthanol,'O2','CO2','H2O','N2'});                      % Find indexes of these species
SpSEthanol=Sp(iSpEthanol);                                                                % Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSpEthanol = length(SpSEthanol);
MiEthanol = [SpSEthanol.Mass];
YfuelEthanol = [1 0 0 0 0];
MEthanol = YfuelEthanol*MiEthanol';
Ethanol_comp = SpSEthanol.Elcomp; %Answer in order [O H C N AR]
% Result gives: [1 6 2 0 0]

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

% C7.76H13.1 + 11.035 O2 −-> 7.76 CO2 + 6.55 H2O chemical reaction for
% Gasoline

% C2H5OH + 3O2 --> 2CO2 + 3 H20 chemical reaction for Ethanol

%% Calculation AF ratio pure gasoline
MolesGasoline = 1/MGasoline; % Amount of moles gasoline per mass unit
MolesOxygenG = MolesGasoline*11.035; % Amount of moles Oxygen for this ... 
% amount of gasoline
MolesNitrogenG = (MolesOxygenG*0.79)/0.21; %Amount of moles Nitrogen for ...
% this amount of gasoline
MassOxygenG = MolesOxygenG*MiGasoline(2); % The mass of this oxygen
MassNitrogenG = MolesNitrogenG*MiGasoline(5);
MassAirG = MassOxygenG + MassNitrogenG; % Mass of air (Oxygen + ...
% Nitrogen)
AirFuelRatioGasoline = MassAirG / 1; % Air−fuel ratio for Gasoline

%% Calculation AF ratio E5
E5_value_Ethanol = 0.05; %volume percentage Ethanol
E5_value_Gasoline = 0.95; %volume percentage Gasoline
%Assume 1 kg of fuel: m_ethanol + m_gasoline = 1

MassGasolineE5 = 1/(1+(E5_value_Ethanol/E5_value_Gasoline)*(DensityEthanol/DensityGasoline));
MassEthanolE5 = 1 - MassGasolineE5;

%Mass air needed for combustion Gasoline in E5
MolesGasolineE5 = MassGasolineE5/MGasoline; %Amount of moles Gasoline in 1kg E5
MolesOxygenE5Gasoline = MolesGasolineE5*11.035; %% Amount of moles Oxygen needed ... 
% for combustion gasoline in E5
MolesNitrogenE5Gasoline = (MolesOxygenE5Gasoline*0.79)/0.21; %% Amount of moles Nitrogen needed ... 
% for combustion gasoline in E5
MassOxygenE5Gasoline = MolesOxygenE5Gasoline*MiGasoline(2); % The mass of this oxygen
MassNitrogenE5Gasoline = MolesNitrogenE5Gasoline*MiGasoline(5); % The mass of this nitrogen

TotalMassAirE5Gasoline = MassOxygenE5Gasoline + MassNitrogenE5Gasoline;

%Mass air needed for combustion Ethanol in E5
MolesEthanolE5 = MassEthanolE5/MEthanol; %Amount of moles Ethanol in 1kg E5
MolesOxygenE5Ethanol = MolesEthanolE5*3; %% Amount of moles Oxygen needed ... 
% for combustion ethanol in E5
MolesNitrogenE5Ethanol = (MolesOxygenE5Ethanol*0.79)/0.21; %% Amount of moles Nitrogen needed ... 
% for combustion ethanol in E5
MassOxygenE5Ethanol = MolesOxygenE5Ethanol*MiEthanol(2);
MassNitrogenE5Ethanol = MolesNitrogenE5Ethanol*MiEthanol(5);

TotalMassAirE5Ethanol = MassOxygenE5Ethanol + MassNitrogenE5Ethanol;

%Total mass air needed for combustion of 1kg E5
TotalMassAirE5 = TotalMassAirE5Gasoline + TotalMassAirE5Ethanol;

%AF ratio for E5
AirFuelRatioE5 = TotalMassAirE5/1;

%% Calculation AF ratio E10
E10_value_Ethanol = 0.10; %volume percentage Ethanol
E10_value_Gasoline = 0.90; %volume percentage Gasoline
%Assume 1 kg of fuel: m_ethanol + m_gasoline = 1

MassGasolineE10 = 1/(1+(E10_value_Ethanol/E10_value_Gasoline)*(DensityEthanol/DensityGasoline));
MassEthanolE10 = 1 - MassGasolineE10;

%Mass air needed for combustion Gasoline in E10
MolesGasolineE10 = MassGasolineE10/MGasoline; %Amount of moles Gasoline in 1kg E10
MolesOxygenE10Gasoline = MolesGasolineE10*11.035; %% Amount of moles Oxygen needed ... 
% for combustion gasoline in E10
MolesNitrogenE10Gasoline = (MolesOxygenE10Gasoline*0.79)/0.21; %% Amount of moles Nitrogen needed ... 
% for combustion gasoline in E10
MassOxygenE10Gasoline = MolesOxygenE10Gasoline*MiGasoline(2); % The mass of this oxygen
MassNitrogenE10Gasoline = MolesNitrogenE10Gasoline*MiGasoline(5); % The mass of this nitrogen

TotalMassAirE10Gasoline = MassOxygenE10Gasoline + MassNitrogenE10Gasoline;

%Mass air needed for combustion Ethanol in E10
MolesEthanolE10 = MassEthanolE10/MEthanol; %Amount of moles Ethanol in 1kg E10
MolesOxygenE10Ethanol = MolesEthanolE10*3; %% Amount of moles Oxygen needed ... 
% for combustion ethanol in E10
MolesNitrogenE10Ethanol = (MolesOxygenE10Ethanol*0.79)/0.21; %% Amount of moles Nitrogen needed ... 
% for combustion ethanol in E10
MassOxygenE10Ethanol = MolesOxygenE10Ethanol*MiEthanol(2);
MassNitrogenE10Ethanol = MolesNitrogenE10Ethanol*MiEthanol(5);

TotalMassAirE10Ethanol = MassOxygenE10Ethanol + MassNitrogenE10Ethanol;

%Total mass air needed for combustion of 1kg E10
TotalMassAirE10 = TotalMassAirE10Gasoline + TotalMassAirE10Ethanol;

%AF ratio for E10
AirFuelRatioE10 = TotalMassAirE10/1;