clear all;
close all;
%% Load the NASA tables

relativepath_to_generalfolder='General'; % relative reference to General folder (assumes the folder is in you working folder)
addpath(relativepath_to_generalfolder); 
% Load Nasadatabase
TdataBase=fullfile('General','Nasa','NasaThermalDatabase');
load(TdataBase);
run("AFcalculations.m");
run("Fraction_calculations.m");
run("system_parameters.m");
run("Carburetor.m")

%% Constants

global Runiv
Runiv=8.314472;

%% Fuel computations

Evalue = 10;          % E-number of the fuel

% Qlvh = Amount of energy per mass of fuel (j)
if Evalue == 0
    Qlhv = 46.4e6;            
elseif Evalue == 5
    Qlhv = 45.58e6;               % Could not find a value on the internet, this is an approximation
elseif Evalue == 10
    Qlhv = 43.54e6;
end



% Composition Ethanol
cFuelEthanol = 'C2H5OH';        %Ethanol
iSpEthanol = myfind({Sp.Name},{cFuelEthanol,'O2','CO2','H2O','N2'});     % Find indexes of these species
SpSEthanol=Sp(iSpEthanol);                                               % Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSpEthanol = length(SpSEthanol);
MiEthanol = [SpSEthanol.Mass];
YfuelEthanol = [1 0 0 0 0];
MEthanol = YfuelEthanol*MiEthanol';         % Molar mass of ethanol

% Composition Air
Xair = [0 0.21 0 0 0.79];                   % Order is important. Note that these are molefractions
MAir = Xair*MiEthanol';                     % Row times Column = inner product 
Yair = Xair.*MiEthanol/MAir;                % Vector. times vector is Matlab's way of making an elementwise multiplication. Note that these are Mass Fractions

%C2H5OH + 3 O2 → 2 CO2 + 3 H2O
MolesEthanol = 1/MEthanol;                  % Amount of moles Ethanol per mass unit
MolesOxygenE = MolesEthanol*3;              % Amount of moles Oxygen for this amount of Ethanol
MassOxygenE = MolesOxygenE*MiEthanol(2);    % The mass of this oxygen
MassAirE = MassOxygenE + (MassOxygenE/Yair(2))*Yair(5); % Mass of air (Oxygen + Nitrogen)
AirFuelRatioEthanol = MassAirE / 1;  % Air-fuel ratio for Ethanol

% Composition Gasoline
cFuelGasoline = 'Gasoline';
iSpGasoline = myfind({Sp.Name},{cFuelGasoline,'O2','CO2','H2O','N2'});                    % Find indexes of these species
SpSGasoline=Sp(iSpGasoline);                                                              % Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSpGasoline = length(SpSGasoline);
MiGasoline = [SpSGasoline.Mass];
YfuelGasoline = [1 0 0 0 0];
MGasoline = YfuelGasoline*MiGasoline';                                                    % Molar mass of gasoline

% 2 C8H18 + 25 O2 → 16 CO2 + 18 H2O     chemical reaction for octane
MolesGasoline = 1/MGasoline;                  % Amount of moles gasoline per mass unit
MolesOxygenG = MolesGasoline*12.5;              % Amount of moles Oxygen for this amount of gasoline
MassOxygenG = MolesOxygenG*MiGasoline(2);    % The mass of this oxygen
MassAirG = MassOxygenG + (MassOxygenG/Yair(2))*Yair(5); % Mass of air (Oxygen + Nitrogen)

VolumeEthanol = Evalue/100;
VolumeGasoline = (100-Evalue)/100;

MassEthanol = DensityEthanol*VolumeEthanol;            % Ethanol = 789 kg/m^3
MassGasoline = DensityGasoline*VolumeGasoline;          % Gasoline = 749 kg/m^3   

%% Initialisation

p(1) = p0;
T(1) = T0;
pad(1) = p(1);
Ca(1) = 1;
V(1) = Vcyl(Ca(1),S,B,l,rc); % Vcyl is a function that computes cylinder volume as function of crank-angle for given B,S,l and rc
m(1) = p(1)*V(1)/Rg_before_comb/T(1);
NCa=720;                % Number of crank-angles
dCa=0.5;                % Stepsize
NSteps=NCa/dCa;
Q_LHV_E0 = LowerHeatingValue(T_ref_QLHV,SpSGasoline,iSpGasoline, MiGasoline);
Q_loss(1) = 0;
Rg = 300;
Cv = 800;

for i = 1:NSteps
    Ca_i = (i - 1) * dCa;  % Current crank angle
    
    % Intake stroke
    if Ca_i >= 0 && Ca_i < 180
        B1(i) = 6.18;
        B2(i) = 0;

    % Compression stroke
    elseif Ca_i >= 180 && Ca_i < 360
        B1(i) = 2.28;
        B2(i) = 0;

    % Power stroke / Combustion / Expansion
    elseif Ca_i >= 360 && Ca_i < 540
        B1(i) = 2.28;
        B2(i) = 3.24 * 10^-3;
        dQcom(i) = massflowFuel * Q_LHV_E0;

    % Exhaust stroke
    elseif Ca_i >= 540 && Ca_i <= 720
        B1(i) = 6.18;
        B2(i) = 0;
    end
end

%% Loop over the crank angles using a For-loop


for i=2:NSteps
        Ca(i)=Ca(i-1)+dCa;
        V(i)=Vcyl(Ca(i),S,B,l,rc); % New volume for current crank-angle
        m(i)=m(i-1); % Mass is constant, valves are closed
        dV=V(i)-V(i-1); % Volume change

        
        dT=(dQcom(i)-p(i-1)*dV)/Cv/m(i-1); % 1st Law dU=dQ-pdV (closed system), adiabatic closed system with constant gas composition and constant Cv
        T(i)=T(i-1)+dT;
        p(i)=m(i)*Rg*T(i)/V(i); % Gaslaw

        % for j=1:NSteps
        %     Cvi(j) = CvNasa(T(i),SpS(j)); % Replaced YourNasa function
        % end
        % Y(i,:) = Y(i-1,:); % Assumes constant mass ratios (needs to be improved)
        % Cv(i) = Y(i,:) * Cvi'; % heat cap at current T

    % Heat loss
    % A(i) = (pi/2)*B^2 + pi*B*(r*cosd(Ca(i)) + sqrt(l^2 - r^2*(sind(Ca(i))^2))); % [m^2] Instantaneous inner cylinder area 
    % p_motor2(i) = (P_ref * (V_ref/V(i))^gamma_comb_out); % [Pa] Motorized cylinder pressure
    % w(i) = B1(i)*S_p + B2(i)*((max(V)*T_ref)/(P_ref*V_ref)) * (p(i) - p_motor2(i)); % [m/s] Effective gas velocity
    % h_woschni(i) = 3.26 * B^(-0.2) * p(i)^(0.8) * T(i)^(-0.55) * w(i)^0.8; % [W/(m^2*K)]
    % Q_loss(i) = h_woschni(i) * A(i) * (T(i) - 330); % [W] Convective heat loss to the inner cylinder wall

end

%% Plot pV-diagram

figure;
plot(Ca, V);
xlabel('Crank angle');
ylabel('Volume (m^3)');
title('Crank angle VC Volume');
grid on;


figure;
plot(V, p);
xlabel('Volume (m^3)');
ylabel('Pressure (Pa)');
title('pV-diagram for the Otto cycle');
grid on;

figure;
loglog(V, p);
xlabel('Volume (m^3)');
ylabel('Pressure (Pa)');
title('pV-diagram for the Otto cycle (Log-Log scale)'); 
grid on;

figure;
plot(Ca, T);
xlabel('Crank angle (Ca)');
ylabel('Temperature(K)');
title('Crank angle over Temperature');
grid on;

% figure;
% plot(Ca, p_motor2);
% xlabel('Crank angle (Ca)');
% ylabel('transfer coefficient h');
% title('Convective heat coefficient vs crank angle (hohenberg)');
% grid on;

%% Function of V_cyl
function V_cyl = Vcyl(Ca, S, B, l, rc)

    % Kinematic relations
    V_d = (pi/4) * B^2 * S; % [m^3] Displacement volume
    V_c = V_d/(rc - 1); % [m^3] Clearance volume
    r  = 1/2 * S; % [m] Lenght of crankshaft

    x = r * cosd(Ca) + sqrt(l^2 - r^2 * sind(Ca).^2); % [m] x postion of piston as a function of theta
    d = l + r - x; % [m] Distance of pistion from TDC as a function of theta

    V_cyl = pi * (B/2)^2 * d + V_c; % [m^3] Free cylinder volume as a function of theta
    
end

%% Function of low heating value
function [Q_LHV] = LowerHeatingValue(T_ref_QLHV,SpSGasoline,iSpGasoline, MiGasoline)

    MoleH2O = 6.55;
    MoleCO2 = 7.76;
    MoleO2 = 11.035;
    MoleN2 = 41.5;
    MassH2O = (MoleH2O*0.0180)/(MiGasoline(1));
    MassCO2 = (MoleCO2*0.0440)/(MiGasoline(1));
    MassO2 = (MoleO2*0.0320)/(MiGasoline(1));
    MassN2 = (MoleN2*0.0280)/(MiGasoline(1));

Q_LHV = -MassN2*HNasa(T_ref_QLHV,SpSGasoline(5))-MassCO2*HNasa(T_ref_QLHV,SpSGasoline(3)) - MassH2O*HNasa(T_ref_QLHV,SpSGasoline(4)) + MassO2*HNasa(T_ref_QLHV,SpSGasoline(2)) + HNasa(T_ref_QLHV,SpSGasoline(1));
end