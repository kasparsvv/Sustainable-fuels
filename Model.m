%% Load the NASA tables

relativepath_to_generalfolder='General'; % relative reference to General folder (assumes the folder is in you working folder)
addpath(relativepath_to_generalfolder); 
% Load Nasadatabase
TdataBase=fullfile('General','Nasa','NasaThermalDatabase');
load(TdataBase);

%% Constants

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
AirFuelRatioGasoline = MassAirG / 1;         % Air-fuel ratio for Gasoline 



VolumeEthanol = Evalue/100;
VolumeGasoline = (100-Evalue)/100;

MassEthanol = 789*VolumeEthanol;            % Ethanol = 789 kg/m^3
MassGasoline = 749*VolumeGasoline;          % Gasoline = 749 kg/m^3   
%Mass fractions
MassFractionEthanol = MassEthanol/(MassEthanol+MassGasoline);
MassFractionGasoline = MassGasoline/(MassEthanol+MassGasoline);

MFuel = MassFractionEthanol*MEthanol + MassFractionGasoline*MGasoline;      % Molar mass of the fuel mixture
mfuel = 1;               %Fuel mass inside the cylinder, used for Qmodel, not defined yet

AirFuelRatio = MassFractionEthanol*AirFuelRatioEthanol + MassFractionGasoline*AirFuelRatioGasoline;     % Air-fuel ratio for the fuel mixture 
MFuelAir = MFuel + AirFuelRatio*MFuel;          % Mass of air = Air-Fuel ratio * Mass of the fuel

Rg = Runiv/MFuelAir;   %Specific gas constant
%Rg = 290;


%% Initialisation

p(1) = p0;
T(1) = T0;
pad(1) = p(1);
Ca(1) = 0;
V(1) = Vcyl(Ca(1),S,B,l,rc); % Vcyl is a function that computes cylinder volume as function of crank-angle for given B,S,l and rc
m(1) = p(1)*V(1)/Rg/T(1);

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
Cv = 1000;

gamma = 1.4;
C1 = p(1)*V(1)^gamma;       % Poisson relations
C2 = T(1)*V(1)^(gamma-1);   % Poisson relations

%% Loop over the crank angles using a For-loop

for i=2:NSteps                          % Calculate values for 1 cycle
    Ca(i)=Ca(i-1)+dCa;
    V(i)=Vcyl(Ca(i),S,B,l,rc);          % New volume for current crank-angle
    m(i)=m(i-1);                        % Mass is constant, valves are closed
    dV=V(i)-V(i-1);                     % Volume change

    % Intake
    if Ca(i) >= 0 && Ca(i) < 180
        p(i) = p0;
    end

    % Compression
    if Ca(i) >= 180 && Ca(i) < 360
        p(i) = C1/V(i)^(gamma);         % Poisson relations
        T(i) = C2/V(i)^(gamma-1);       % Poisson relations
    end

    % Ignition
    if Ca(i) == 360
        dQcom(i) = 730;                         % Heat Release by combustion
        dT(i)=(dQcom(i)-p(i-1)*dV)/Cv/m(i-1);   % 1st Law dU=dQ-pdV (closed system)
        T(i)=T(i-1)+dT(i);
        p(i)=m(i)*Rg*T(i)/V(i);                 % Gaslaw

        % for n=1:NSp
        % Cvi(n) =CvNasa(T(i),SpS(n));           % Get Cv from Nasa-table
        % end
        % Cv = Yi_after * Cvi';

        % for n=1:NSp
        % Cpi(n) =CpNasa(T(i),SpS(n));           % Get Cp from Nasa-table
        % end
        % Cp = Yi_after * Cvi';

        % gamma = Cp/Cv;

        C3 = p(i)*V(i)^gamma;                   % Poisson relations
        C4 = T(i)*V(i)^(gamma-1);               % Poisson relations
    end

    % Power stroke
    if Ca(i) > 360 && Ca(i) < 540
        p(i) = C3/V(i)^(gamma);         % Poisson relations
        T(i) = C4/V(i)^(gamma-1);       % Poisson relations
    end

    % Heat release
    if Ca(i) == 540      

        m(i) = m(i-1);     

        p(i) = p0;         
        T(i) = T0;  
        Q_c = Cv * m(i) * (T(i)-T(i-1));
    end

    % Exhaust
    if Ca(i) >= 540 && Ca(i) <= 720
        p(i) = p0;
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