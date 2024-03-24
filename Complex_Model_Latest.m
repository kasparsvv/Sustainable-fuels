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


%% Determine load and E-value

Evalue = 0;          % E-number of the fuel
Load = 0;            %Determine load, either full (1), half (0.5) or no (0) load

% Qlvh = Amount of energy per mass of fuel (j)
if Evalue == 0
    Qlhv = 46.4e6;
    if Load == 1
        p0 = 101235;            
    elseif Load== 0.5
        p0 = 61300;               
    elseif Load == 0
        p0 = 27000;
    end
elseif Evalue == 5
    Qlhv = 45.58e6;               % Could not find a value on the internet, this is an approximation
    if Load == 1
        p0 = 101235;            
    elseif Load== 0.5
        p0 = 61200;               
    elseif Load == 0
        p0 = 27000;
    end
elseif Evalue == 10
    Qlhv = 43.54e6;
    if Load == 1
        p0 = 101235;            
    elseif Load== 0.5
        p0 = 61155;               
    elseif Load == 0
        p0 = 26912;
    end
end





%% Fuel computations


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

%% Initialisation

p(1) = p0;
T(1) = T0;
pad(1) = p(1);
Ca(1) = 0;
V(1) = Vcyl(Ca(1),S,B,l,rc); % Vcyl is a function that computes cylinder volume as function of crank-angle for given B,S,l and rc
m(1) = p(1)*V(1)/Rg_before_comb/T(1);
Q_loss(1) = 0;
gamma_comb_out = 1.236322380761674;
gamma_comb_in = 1.360214357127973;

Rg_before_comb = 2.744179546614579e+02;
Rg_after_comb = 2.860747563192256e+02;


for i = 1:NSteps+1
    Ca_i = (i - 1) * dCa;  % Current crank angle

    % Intake stroke
    if Ca_i >= 0 && Ca_i < 180
        B1(i) = 6.18;
        B2(i) = 0;
        
    % Compression stroke
    elseif Ca_i >= 180 && Ca_i < 350
        B1(i) = 2.28;
        B2(i) = 0;
        
    % Power stroke / Combustion / Expansion
    elseif Ca_i >= 350 && Ca_i < 540
        B1(i) = 2.28;
        B2(i) = 3.24 * 10^-3;
        
    % Exhaust stroke
    elseif Ca_i >= 540 && Ca_i <= 720
        B1(i) = 6.18;
        B2(i) = 0;
    end
end


%% Reference state + intake
for i=2:NSteps+1
    Ca(i)=Ca(i-1)+dCa;
    V(i)=Vcyl(Ca(i),S,B,l,rc); 
    dV=V(i)-V(i-1); 

    % Intake
    if Ca(i) >= 0 && Ca(i) < 180
        p(i) = p0;
        T(i) = T0;
        m(i) = p(i)*V(i)/Rg_before_comb/T(i);
    end
       
    % Reference state
    if Ca(i) == 180
        P_ref = p(360);
        T_ref = T(360);
        V_ref = V(360);
    end
end

%% Loop over the crank angles using a For-loop
for i=2:NSteps+1
    Ca(i)=Ca(i-1)+dCa;
    V(i)=Vcyl(Ca(i),S,B,l,rc); 
    dV=V(i)-V(i-1); 

    for n=1:5
            Cvi(n) = CvNasa(T(360),SpSGasoline(n));
            Cpi(n) = CpNasa(T(360),SpSGasoline(n));
    end
    Cv_comp_in = dot(Y_comp_in,Cvi);
    Cp_comp_in = dot(Y_comp_in,Cpi);


    % Compression
    if Ca(i) >= 180 && Ca(i) < 350
        m(i) = p(1)*V(361)/(Rg_before_comb*T(1));
        dT(i)=(-Q_loss(i-1) -p(i-1)*dV)/Cv_comp_in/m(i-1);
        T(i)=T(i-1)+dT(i);
        p(i)=m(i)*Rg_before_comb*T(i)/V(i);  

    end

    % Ignition
    if Ca(i) >= 350 && Ca(i) <= 540

        for n=1:5
        Cvi_comb_in(n) =CvNasa(T(720),SpSGasoline(n));
        end
        Cv_comb_in = dot(Y_comp_in,Cvi_comb_in);
        m(i) = m(i-1);        
        m_fuel = m(365)/(1+AirFuelRatioGasoline);
        Q_LHV_E0 = LowerHeatingValue(T_ref_QLHV,SpSGasoline,iSpGasoline, MiGasoline);
        dQcomb(i) = QModel(Ca(i),CaS,CaD,m_fuel,Q_LHV_E0);

        dT(i)=(dQcomb(i) - Q_loss(i-1) - p(i-1)*dV)/Cv_comb_in/m(i);
        T(i)=T(i-1)+dT(i);
        p(i)=m(i)*Rg_before_comb*T(i)/V(i); 
      

    for n=1:5
        Cvi_comb_out(n) = CvNasa(T(721),SpSGasoline(n));
        Cpi_comb_out(n) = CpNasa(T(721),SpSGasoline(n));
    end
    Cv_comb_out = dot(Y_comb_out,Cvi_comb_out);
    Cp_comb_out = dot(Y_comb_out,Cpi_comb_out);

    end

    for n=1:5
        Cvi_ps_out(n) =CvNasa(T(1080),SpSGasoline(n));
    end
    Cv_ps_out = dot(Y_comb_out,Cvi_ps_out);

    % Heat release
    if Ca(i) == 540      
        m(i) = p(1)*V(361)/(Rg_after_comb*T(1));    
        p(i) = p(i-1);
        T(i) = T(i-1);
        Q_c = Cv_ps_out * m(i) * (T(i)-T(i-1));
    end

    % Exhaust
    if Ca(i) >= 540 && Ca(i) <= 720
        p(i) = p(i-1);
        T(i) = T(i-1);
        m(i) = p(i)*V(i)/Rg_after_comb/T(i);
    end

    % Closing cycle
    if Ca(i) >= 720
        p(i) = p0;
        T(i) = T0;
        m(i) = p(i)*V(i)/Rg_before_comb/T(i);
    end

    % Heat loss
    A(i) = (pi/2)*B^2 + pi*B*(r*cosd(Ca(i)) + sqrt(l^2 - r^2*(sind(Ca(i))^2))); % [m^2] Instantaneous inner cylinder area 
    p_motor2(i) = (P_ref * (V_ref/V(i))^gamma_comb_out); % [Pa] Motorized cylinder pressure
    w(i) = B1(i)*S_p + B2(i)*(V_d*T_ref)/(P_ref*V_ref) * (p(i) - p_motor2(i)); % [m/s] Effective gas velocity
    h_woschni(i) = 3.26 * B^(-0.2) * (p(i)/1000)^(0.8) * T(i)^(-0.55) * w(i)^0.8; % [W/(m^2*K)]
    Q_loss(i) = h_woschni(i) * A(i) * (T(i) - 450); % [W] Convective heat loss to the inner cylinder wall
    Q_loss(i) = Q_loss(i)/360/50; % [W] Convective heat loss to the inner cylinder wall

end


heat = sum(dQcomb);
%% Efficiency and Power Calculations
RPM = 3000; % rounds per minute
n_rpc = 2; %number of rounds per cycle

W_E0= trapz(V,p);
totaldQcomb = sum(dQcomb);
efficiency = W_E0/totaldQcomb*100;
P_E0 = W_E0 * (RPM/60) * (1/n_rpc);
bsfc = m_fuel*1000/(W_E0/3600000);

%% Plot pV-diagram

% figure;
% plot(Ca, V);
% xlabel('Crank angle');
% ylabel('Volume (m^3)');
% title('Crank angle VC Volume');
% grid on;
% 
figure;
plot(Ca, Q_loss);
xlabel('Crank angle');
ylabel('Volume (m^3)');
title('heat loss');
grid on;
% 
% figure;
% plot(Ca, p_motor2);
% xlabel('Crank angle');
% ylabel('Volume (m^3)');
% title('heat loss');
% grid on;


figure;
plot(V, p);
xlabel('Volume (m^3)');
ylabel('Pressure (Pa)');
title('pV-diagram for the complex cycle');
grid on;

% figure;
% loglog(V, p);
% xlabel('Volume (m^3)');
% ylabel('Pressure (Pa)');
% title('pV-diagram for the complex cycle (Log-Log scale)'); 
% grid on;
% 
figure;
plot(Ca, T);
xlabel('Crank angle (Ca)');
ylabel('Temperature(K)');
xlim([180; 540])
title('Crank angle over Temperature');
grid on;

figure;
plot(Ca, h_woschni);
xlabel('Crank angle (Ca)');
ylabel('transfer coefficient h');
xlim([180; 540])
title('Convective heat coefficient vs crank angle (WOSCHNI)');
grid on;

% figure;
% plot(Ca, h_hohenberg);
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

%% Wiebe function

function [dQcomb] = QModel(Ca,CaS,CaD,mfuel,Q_LHV_E0)

    global Runiv
    if (isempty(Runiv))
        fprintf('[Qmodel] Assign global Runiv\n');
        return
    end

    a = 5; % Wiebe constant 
    n = 3; % Wiebe constant

    xb = 1 - exp(-a * ((Ca - CaS) / CaD) .^ n); % Fuel consumption based on the crank angle
    dQcomb = Q_LHV_E0 * mfuel * n * a * ((1 - xb) / CaD) .* ((Ca - CaS) / CaD) .^ (n - 1); % Heat release, Formula from project handbook page 14

end