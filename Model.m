clear all;
%% Load the NASA tables

relativepath_to_generalfolder='General'; % relative reference to General folder (assumes the folder is in you working folder)
addpath(relativepath_to_generalfolder); 
% Load Nasadatabase
TdataBase=fullfile('General','Nasa','NasaThermalDatabase');
load(TdataBase);
run("AFcalculations.m");
run("Fraction_calculations.m");
run("system_parameters.m")

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

P_ref = 1.859606618881350e+06; % [Pa] Reference pressure taken at crangle angle of 360
V_ref = 2.640285571312627e-05; % [Pa] Reference volume  taken at crangle angle of 360
T_ref = 6.333244594935253e+02; % [Pa] Reference temperature taken at crangle angle of 360

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


for i=2:NSteps                          % Calculate values for 1 cycle
    Ca(i)=Ca(i-1)+dCa;
    V(i)=Vcyl(Ca(i),S,B,l,rc);          % New volume for current crank-angle
    dV=V(i)-V(i-1);                     % Volume change

    % Intake
    if Ca(i) >= 0 && Ca(i) < 180
        p(i) = p0;
        T(i) = T0;
        m(i) = p(i)*V(i)/Rg_before_comb/T(i);

    end
    for n=1:5
            Cvi(n) = CvNasa(T(360),SpSGasoline(n));
            Cpi(n) = CpNasa(T(360),SpSGasoline(n));
    end
    %Cv_comp_in = Y_comp_in.*Cvi'
    %Cp_comp_in = Y_comp_in.*Cpi'

    Cv_comp_in = dot(Y_comp_in,Cvi);
    Cp_comp_in = dot(Y_comp_in,Cpi);

    gamma_comp_in = Cp_comp_in/Cv_comp_in;


    % Compression
    if Ca(i) >= 180 && Ca(i) < 360
        C1 = p(360)*V(360)^gamma_comp_in;
        C2 = T(360)*V(360)^(gamma_comp_in-1);

        m(i) = p(1)*V(361)/(Rg_before_comb*T(1));
        p(i) = C1/V(i)^(gamma_comp_in);         % Poisson relations
        T(i) = C2/V(i)^(gamma_comp_in - 1);       % Poisson relations
    end

    % Ignition
    if Ca(i) == 360
        for n=1:5
        Cvi_comb_in(n) =CvNasa(T(720),SpSGasoline(n));           % Get Cv from Nasa-table
        end
        Cv_comb_in = dot(Y_comp_in,Cvi_comb_in);
        m(i) = p(1)*V(361)/(Rg_before_comb*T(1));
        dQcom(i) = 730;                         % Heat Release by combustion
        dT(i)=(dQcom(i)-p(i-1)*dV)/Cv_comb_in/m(i-1);   % 1st Law dU=dQ-pdV (closed system)
        T(i)=T(i-1)+dT(i);
        p(i)=m(i)*Rg_before_comb*T(i)/V(i);                 % Gaslaw

        
        % for n=1:NSp
        % Cpi(n) =CpNasa(T(i),SpS(n));           % Get Cp from Nasa-table
        % end
        % Cp = Yi_after * Cvi';

        % gamma = Cp/Cv;

        
    end
    for n=1:5
    Cvi_comb_out(n) = CvNasa(T(721),SpSGasoline(n));
    Cpi_comb_out(n) = CpNasa(T(721),SpSGasoline(n));
    end
    %Cv_comp_in = Y_comp_in.*Cvi'
    %Cp_comp_in = Y_comp_in.*Cpi'

    Cv_comb_out = dot(Y_comb_out,Cvi_comb_out);
    Cp_comb_out = dot(Y_comb_out,Cpi_comb_out);

    gamma_comb_out = Cp_comb_out/Cv_comb_out;

    % Power stroke
    if Ca(i) > 360 && Ca(i) < 540
        m(i) = p(1)*V(361)/(Rg_before_comb*T(1));
        C3 = p(721)*V(721)^gamma_comb_out;
        C4 = T(721)*V(721)^(gamma_comb_out-1);
        p(i) = C3/V(i)^(gamma_comb_out);         % Poisson relations
        T(i) = C4/V(i)^(gamma_comb_out-1);       % Poisson relations
    end
    for n=1:5
        Cvi_ps_out(n) =CvNasa(T(1080),SpSGasoline(n));           % Get Cv from Nasa-table
    end
    Cv_ps_out = dot(Y_comb_out,Cvi_ps_out);
    % Heat release
    if Ca(i) == 540      

        m(i) = p(1)*V(361)/(Rg_before_comb*T(1));     

        p(i) = p0;         
        T(i) = T0;  
        Q_c = Cv_ps_out * m(i) * (T(i)-T(i-1));
    end

    % Exhaust
    if Ca(i) >= 540 && Ca(i) <= 720
        p(i) = p0;
        T(i) = T0;
        m(i) = p(i)*V(i)/Rg_after_comb/T(i);
    end


    % A(i) = (pi/2)*B^2 + pi*B*(r*cosd(Ca(i)) + sqrt(l^2 - r^2*(sind(Ca(i))^2)));
    % Q_loss(i) = h_woschni(i) * A(i) * (T(i) - 330);

    % p_motor1(i) = (((rc *  max(V)/(rc - 1))^gamma_comb_out * P_atm)/V(i)^gamma_comb_out); % [Pa] motorized cylinder pressure

    p_motor2(i) = (P_ref * (V_ref/V(i))^gamma_comb_out);

    w(i) = B1(i)*S_p + B2(i)*((max(V)*T_ref)/(P_ref*V_ref)) * (p(i) - p_motor2(i)); % [m/s] Effective gas velocity

    h_woschni(i) = 3.26 * B^(-0.2) * p(i)^(0.8) * T(i)^(-0.55) * w(i)^0.8; % [W/(m^2*K)]
    h_hohenberg(i) = 140 * V(i)^(-0.06) * p(i)^(0.8) * T(i)^(-0.4) * (S_p + 1.4)^(0.8);
    h_eichelberg(i) = 7.799 * 10^(-3) * S_p^(1/3) * p(i)^(0.5) * T(i)^(0.5);

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

figure;
plot(Ca, h_woschni);
xlabel('Crank angle (Ca)');
ylabel('transfer coefficient h');
title('Convective heat coefficient vs crank angle (WOSCHNI)');
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