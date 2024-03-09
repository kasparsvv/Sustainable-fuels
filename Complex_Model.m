
%% Load the NASA tables

relativepath_to_generalfolder='General'; % relative reference to General folder (assumes the folder is in you working folder)
addpath(relativepath_to_generalfolder); 
% Load Nasadatabase
TdataBase=fullfile('General','Nasa','NasaThermalDatabase');
load(TdataBase);
run("AFcalculations.m");
run("Fraction_calculations.m");
run("system_parameters.m")
iSp = myfind({Sp.Name},{'Gasoline','O2','CO2','H2O','N2'});                      % Find indexes of these species
SpS = Sp(iSp);                                                                % Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSp = length(SpS);

%% Constants

global Runiv
Runiv=8.314472;

gamma_comb_out = 1.4;

NCa=720; % Number of crank-angles
dCa=0.5; % Stepsize
NSteps=NCa/dCa;

mfuel = 10;
Cv = 1000;

CaS = -15; % Crank angle of ignition timing
CaD = 35; % Crank angle combustion duration

Q_LHV_g = 43.2e6; %LHV of gasoline
Q_LHV_e = 26.8e6; %LHV of ethanol

% Qcom_calc = Q_LHV * m_fuel;
Qcom = 0;

S_p = 2 * RPM/60 * S;   % [m/s] Mean piston speed

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

%% Initialisation

p(1) = p0;
T(1) = T0;
pad(1) = p(1);
Ca(1) = 0;
V(1) = Vcyl(Ca(1),S,B,l,rc); % Vcyl is a function that computes cylinder volume as function of crank-angle for given B,S,l and rc
m(1) = p(1)*V(1)/Rg_before_comb/T(1);
Q_loss(1) = 0;
Rg = 10;
dQcomb(1) = 0;

%% Constant arrays

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
        
    % Exhaust stroke
    elseif Ca_i >= 540 && Ca_i <= 720
        B1(i) = 6.18;
        B2(i) = 0;
    end
end
%% Loop over crank-angle, with 'for' construction


for i=2:NSteps

    Ca(i)=Ca(i-1)+dCa;
    V(i)=Vcyl(Ca(i),S,B,l,rc); % New volume for current crank-angle
    m(i)=m(i-1); % Mass is constant, valves are closed
    dV=V(i)-V(i-1); % Volume change

    dQcom(i) = QModel(Ca(i),CaS,CaD,mfuel,Q_LHV_g); % Heat Release by combustion
    Qcom = Qcom + dQcom(i-1);

    dQ = dQcom(i);
    dT=(dQ-p(i-1)*dV)/Cv/m(i-1); % 1st Law dU=dQ-pdV (closed system)
    T(i)=T(i-1)+dT;
    p(i)=m(i)*Rg_before_comb*T(i)/V(i); % Gaslaw

    A(i) = (pi/2)*B^2 + pi*B*(r*cosd(Ca(i)) + sqrt(l^2 - r^2*(sind(Ca(i))^2))); % [m^2] Instantaneous inner cylinder area 
    p_motor2(i) = (P_ref * (V_ref/V(i))^gamma_comb_out); % [Pa] Motorized cylinder pressure
    w(i) = B1(i)*S_p + B2(i)*((max(V)*T_ref)/(P_ref*V_ref)) * (p(i) - p_motor2(i)); % [m/s] Effective gas velocity
    h_woschni(i) = 3.26 * B^(-0.2) * p(i)^(0.8) * T(i)^(-0.55) * w(i)^0.8; % [W/(m^2*K) Convective heat coefficient
    Q_loss(i) = h_woschni(i) * A(i) * (T(i) - 200); % [W] Convective heat loss to the inner cylinder wall

        % Intake
    if Ca(i) >= 0 && Ca(i) <= 180
        p(i) = p0;
        T(i) = T0;
        m(i) = p(i)*V(i)/Rg_before_comb/T(i);
    end

        % Exhaust
    if Ca(i) >= 540 && Ca(i) <= 720
        p(i) = p0;
        T(i) = T0;
        m(i) = p(i)*V(i)/Rg_after_comb/T(i);
    end


end

%% Plot pV-diagram

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

figure;
plot(Ca, h_woschni);
xlabel('Crank angle (Ca)');
ylabel('transfer coefficient h');
title('Convective heat coefficient vs crank angle (WOSCHNI)');
grid on;

figure;
plot(Ca, Q_loss);
xlabel('Crank angle (Ca)');
ylabel('Heat loss [W]');
title('Heat loss vs crank angle');
grid on;

figure;
plot(Ca, h_hohenberg);
xlabel('Crank angle (Ca)');
ylabel('transfer coefficient h');
title('Convective heat coefficient vs crank angle (hohenberg)');
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

function [Q_LHV] = LowerHeatingValue(T_ref_QLHV, SpSGasoline, iSpGasoline, MiGasoline)

    MoleH2O = 6.55;
    MoleCO2 = 7.76;
    MoleO2 = 11.035;
    MoleN2 = 41.5;
    MassH2O = (MoleH2O*0.0180)/(MiGasoline(1));
    MassCO2 = (MoleCO2*0.0440)/(MiGasoline(1));
    MassO2 = (MoleO2*0.0320)/(MiGasoline(1));
    MassN2 = (MoleN2*0.0280)/(MiGasoline(1));

    Q_LHV = -MassN2*HNasa(T_ref_QLHV, SpSGasoline(5)) - MassCO2*HNasa(T_ref_QLHV, SpSGasoline(3)) - MassH2O*HNasa(T_ref_QLHV, SpSGasoline(4)) + MassO2*HNasa(T_ref_QLHV, SpSGasoline(2)) + HNasa(T_ref_QLHV, SpSGasoline(1));

    % Call QModel function passing Q_LHV
    QModel(Q_LHV);
end

function [dQcomb] = QModel(Ca,CaS,CaD,mfuel,Q_LHV)

    global Runiv
    if (isempty(Runiv))
        fprintf('[Qmodel] Assign global Runiv\n');
        return
    end

    a = 5; % Wiebe constant 
    n = 3; % Wiebe constant

    xb = 1 - exp(-a * ((Ca - CaS) / CaD) .^ n); % Fuel consumption based on the crank angle
    dQcomb = Q_LHV * mfuel * n * a * ((1 - xb) / CaD) .* ((Ca - CaS) / CaD) .^ (n - 1); % Heat release, Formula from project handbook page 14

end
