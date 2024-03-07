clear all;
%% Load the NASA tables

relativepath_to_generalfolder='General'; % relative reference to General folder (assumes the folder is in you working folder)
addpath(relativepath_to_generalfolder); 
% Load Nasadatabase
TdataBase=fullfile('General','Nasa','NasaThermalDatabase');
load(TdataBase);
run("AFcalculations.m");
run("Fraction_calculations.m");

%% Constants

global Runiv
Runiv=8.314472;

% Crank angle steps
NCa=720;                % Number of crank-angles
dCa=0.5;                % Stepsize
NSteps=NCa/dCa;

% Create arrays
p = zeros(1, NSteps);
V = zeros(1, NSteps);
T = zeros(1, NSteps);
dQcom = zeros(1, NSteps);
dT = zeros(1, NSteps);
m = zeros(1, NSteps);

% System constant
p0 = 1.01235e5 ;        % Ambient Pressure (pa)
T0 = 293;               % Ambient Temparature (K)
rc = 8.5;               % Compression ratio (-)
Runiv = 8.314472;       % Universal gas constant (J / mol·K)
CaS = 210;              % Crank angle at start of combustion
CaD = 225-210;          % Combustion duration (in terms of crank angle)
n = 3;                  % Wiebe form factor, Project handbook says 3 is often used
a = 5;                  % Wiebe efficiency factor, Project handbook says 5 is often used

S = 0.055;              % [m] Stroke
r  = 1/2 * S;           % [m] Length of crankshaft
l = 0.0842;             % [m] Length of the connecting rod
B = 0.0677;             % [m] Bore dimension
TDC = 0.003;            % [m] Top dead center
BDC = S + TDC;          % [m] Bottom dead center
V_c = pi/4 * B^2 * TDC; % [m63] Compression volume (will change)
P_atm  = 1;             % [Bar] Atmospheric pressure (assumed for now)
Ca(1) = 0;              % Initial crank angle

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
%AirFuelRatioGasoline = MassAirG / 1;         % Air-fuel ratio for Gasoline 



VolumeEthanol = Evalue/100;
VolumeGasoline = (100-Evalue)/100;

MassEthanol = DensityEthanol*VolumeEthanol;            % Ethanol = 789 kg/m^3
MassGasoline = DensityGasoline*VolumeGasoline;          % Gasoline = 749 kg/m^3   
%Mass fractions
%MassFractionEthanol = MassEthanol/(MassEthanol+MassGasoline);
%MassFractionGasoline = MassGasoline/(MassEthanol+MassGasoline);

%MFuel = MassFractionEthanol*MEthanol + MassFractionGasoline*MGasoline;      % Molar mass of the fuel mixture
%mfuel = 1;               %Fuel mass inside the cylinder, used for Qmodel, not defined yet

%AirFuelRatio = MassFractionEthanol*AirFuelRatioEthanol + MassFractionGasoline*AirFuelRatioGasoline;     % Air-fuel ratio for the fuel mixture 
%MFuelAir = MFuel + AirFuelRatio*MFuel;          % Mass of air = Air-Fuel ratio * Mass of the fuel

%Rg = Runiv/MFuelAir;   %Specific gas constant
%Rg = 290;
%% Reference temperature for LHV
T_ref_QLHV = 20+273.15;

%% Initialisation

p(1) = p0;
T(1) = T0;
pad(1) = p(1);
Ca(1) = 0;
V(1) = Vcyl(Ca(1),S,B,l,rc); % Vcyl is a function that computes cylinder volume as function of crank-angle for given B,S,l and rc
m(1) = p(1)*V(1)/Rg_E5_before_comb/T(1);

%% Poisson relations

% for n=1:NSp
%     Cvi(n) =CvNasa(T(1),SpS(n));           % Get Cv from Nasa-table
% end
% Cv = Yi_before * Cvi'; 

% for n=1:NSp
%     Cpi(n) =CpNasa(T(1),SpS(n));           % Get Cp from Nasa-table
% end
% Cp = Yi_before * Cpi';

% gamma = Cp/Cv;
%Cv = 1000;

%gamma = 1.4;
%C1 = p(1)*V(1)^gamma;       % Poisson relations
%C2 = T(1)*V(1)^(gamma-1);   % Poisson relations

%% Loop over the crank angles using a For-loop

for i=2:NSteps                          % Calculate values for 1 cycle
    Ca(i)=Ca(i-1)+dCa;
    V(i)=Vcyl(Ca(i),S,B,l,rc);          % New volume for current crank-angle
    dV=V(i)-V(i-1);                     % Volume change

    % Intake
    if Ca(i) >= 0 && Ca(i) < 180
        p(i) = p0;
        T(i) = T0;
        m(i) = p(i)*V(i)/Rg_E5_before_comb/T(i);

    end
    for n=1:6
            Cvi(n) = CvNasa(T(360),SpSE5(n));
            Cpi(n) = CpNasa(T(360),SpSE5(n));
    end
    %Cv_comp_in = Y_comp_in.*Cvi'
    %Cp_comp_in = Y_comp_in.*Cpi'

    Cv_comp_in = dot(Y_E5_comp_in,Cvi);
    Cp_comp_in = dot(Y_E5_comp_in,Cpi);

    gamma_comp_in = Cp_comp_in/Cv_comp_in;


    % Compression
    if Ca(i) >= 180 && Ca(i) < 360
        C1 = p(360)*V(360)^gamma_comp_in;
        C2 = T(360)*V(360)^(gamma_comp_in-1);

        m(i) = p(1)*V(361)/(Rg_E5_before_comb*T(1));
        p(i) = C1/V(i)^(gamma_comp_in);         % Poisson relations
        T(i) = C2/V(i)^(gamma_comp_in - 1);       % Poisson relations
    end

    % Ignition
    if Ca(i) == 360
        for n=1:6
        Cvi_comb_in(n) =CvNasa(T(720),SpSE5(n));           % Get Cv from Nasa-table
        end
        Cv_comb_in = dot(Y_E5_comp_in,Cvi_comb_in);
        Q_LHV_E5 = LowerHeatingValue(T_ref_QLHV,SpSGasoline,SpSEthanol, MassH20E5Gasoline, MassH20E5Ethanol,MassCO2E5Gasoline,MassCO2E5Ethanol, MassOxygenE5Gasoline, MassOxygenE5Ethanol,MassNitrogenE5Gasoline,MassNitrogenE5Ethanol,MassGasolineE5, MassEthanolE5);
        m(i) = p(1)*V(361)/(Rg_E5_before_comb*T(1));
        m_fuel = m(365)/(1+AirFuelRatioE5);
        dQcom = m_fuel*Q_LHV_E5;                         % Heat Release by combustion
        dT(i)=(dQcom-p(i-1)*dV)/Cv_comb_in/m(i-1);   % 1st Law dU=dQ-pdV (closed system)
        T(i)=T(i-1)+dT(i);
        p(i)=m(i)*Rg_E5_before_comb*T(i)/V(i);                 % Gaslaw

        
        % for n=1:NSp
        % Cpi(n) =CpNasa(T(i),SpS(n));           % Get Cp from Nasa-table
        % end
        % Cp = Yi_after * Cvi';

        % gamma = Cp/Cv;

        
    end
    for n=1:6
    Cvi_comb_out(n) = CvNasa(T(721),SpSE5(n));
    Cpi_comb_out(n) = CpNasa(T(721),SpSE5(n));
    end
    %Cv_comp_in = Y_comp_in.*Cvi'
    %Cp_comp_in = Y_comp_in.*Cpi'

    Cv_comb_out = dot(Y_E5_comb_out,Cvi_comb_out);
    Cp_comb_out = dot(Y_E5_comb_out,Cpi_comb_out);

    gamma_comb_out = Cp_comb_out/Cv_comb_out;

    % Power stroke
    if Ca(i) > 360 && Ca(i) < 540
        m(i) = p(1)*V(361)/(Rg_E5_before_comb*T(1));
        C3 = p(721)*V(721)^gamma_comb_out;
        C4 = T(721)*V(721)^(gamma_comb_out-1);
        p(i) = C3/V(i)^(gamma_comb_out);         % Poisson relations
        T(i) = C4/V(i)^(gamma_comb_out-1);       % Poisson relations
    end
    for n=1:6
        Cvi_ps_out(n) =CvNasa(T(1080),SpSE5(n));           % Get Cv from Nasa-table
    end
    Cv_ps_out = dot(Y_E5_comb_out,Cvi_ps_out);
    % Heat release
    if Ca(i) == 540      

        m(i) = p(1)*V(361)/(Rg_E5_before_comb*T(1));     

        p(i) = p0;         
        T(i) = T0;  
        Q_c = Cv_ps_out * m(i) * (T(i)-T(i-1));
    end

    % Exhaust
    if Ca(i) >= 540 && Ca(i) <= 720
        p(i) = p0;
        T(i) = T0;
        m(i) = p(i)*V(i)/Rg_E5_after_comb/T(i);
    end
    
end

%% Efficiency and Power Calculations
RPM = 3000; % rounds per minute
n_rpc = 2; %number of rounds per cycle

W_E5= trapz(V,p);
efficiency = W_E5/dQcom*100;
gamma_average = (gamma_comb_out+gamma_comp_in) / 2;
ottoefficiency =(1-(1/rc)^(gamma_average-1)) *100;
P_E5 = W_E5 * (RPM/60) * (1/n_rpc);
bsfc = m_fuel*1000/(W_E5/3600000);
%% Plot pV-diagram

figure;
plot(V, p);
xlabel('Volume (m^3)');
ylabel('Pressure (Pa)');
title('pV-diagram for the Otto cycle E5');
grid on;

figure;
loglog(V, p);
xlabel('Volume (m^3)');
ylabel('Pressure (Pa)');
title('pV-diagram for the Otto cycle E5 (Log-Log scale)'); 
grid on;


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
function [Q_LHV] = LowerHeatingValue(T_ref_QLHV,SpSGasoline,SpSEthanol, MassH20E5Gasoline, MassH20E5Ethanol,MassCO2E5Gasoline,MassCO2E5Ethanol, MassOxygenE5Gasoline, MassOxygenE5Ethanol,MassNitrogenE5Gasoline,MassNitrogenE5Ethanol,MassGasolineE5, MassEthanolE5)

    MassH2O = MassH20E5Gasoline + MassH20E5Ethanol;
    MassCO2 = MassCO2E5Gasoline + MassCO2E5Ethanol;
    MassO2 = MassOxygenE5Gasoline + MassOxygenE5Ethanol;
    MassN2 = MassNitrogenE5Gasoline + MassNitrogenE5Ethanol;

Q_LHV = -MassCO2*HNasa(T_ref_QLHV,SpSGasoline(3)) - MassH2O*HNasa(T_ref_QLHV,SpSGasoline(4)) + MassGasolineE5*HNasa(T_ref_QLHV,SpSGasoline(1)) + MassEthanolE5*HNasa(T_ref_QLHV,SpSEthanol(1));
end