%% Parameter values. These are retrieved in all other scripts.
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
V_d = (pi/4) * B^2 * S; % [m^3] Displacement volume
rc = 10;                % Compression ratio (assumed for now)
V_c = V_d/(rc - 1);     % [m^3] Clearance volume
P_atm  = 1;             % [Bar] Atmospheric pressure (assumed for now)