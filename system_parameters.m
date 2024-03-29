%% System constants

p0 = 1.01235e5 ;        % Ambient Pressure (pa)
T0 = 293;               % Ambient Temparature (K)
rc = 8.5;               % Compression ratio (-)
Runiv = 8.314472;       % Universal gas constant (J / mol·K)
n = 3;                  % Wiebe form factor, Project handbook says 3 is often used
a = 5;                  % Wiebe efficiency factor, Project handbook says 5 is often used

S = 0.055;              % [m] Stroke
r  = 1/2 * S;           % [m] Length of crankshaft
l = 0.0842;             % [m] Length of the connecting rod
B = 0.0677;             % [m] Bore dimension
TDC = 0.003;            % [m] Top dead center
BDC = S + TDC;          % [m] Bottom dead center
V_c = pi/4 * B^2 * TDC; % [m^3] Compression volume
V_d = pi/4 * B^2  * S;  % [m^3] Swept/Displacement volume
P_atm  = 1;             % [Bar] Atmospheric pressure (assumed for now)
Ca(1) = 0;              % Initial crank angle
RPM = 3000;             % [Hz] Rotation per minute of crankshaft

S_p = 2 * RPM/60 * S;   % [m/s] Mean piston speed
CaD = 60;               % Combustion duration
CaS = 350;              % Crank angle at start of combustion
CaF = CaS + CaD;        % Final crank angle at end of combustion
%% Loop variables/presets

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

%% Reference temperature for LHV
T_ref_QLHV = 20+273.15;