%% Defining constants

S = 0.058; % [m] Stroke
r  = 1/2 * S; % [m] Lenght of crankshaft
l = 0.0842; % [m] Lenght of the connecting rod
B = 0.0677; % [m] Bore dimension
TDC = 0.003; % [m] Top dead center
BDC = S + TDC; % [m] Bottom dead center
V_c = pi/4 * B^2 * TDC; % [m] Compression volume
P_atm  = 1; % [Bar] Atmospheric pressure

%% Defining the kinematic equations of the engine as functions of theta

theta_deg = linspace(-15, 215, 231); % Intake stroke in degrees. Lasts from -15 to 215
theta_rad = deg2rad(theta_deg); % Intake stroke in radians for ease of use

x = r * cos(theta_rad) + sqrt(l^2 - r^2 * sin(theta_rad).^2); % [m] x postion of piston as a function of theta
d = l + r - x; % [m] Distance of pistion from TDC as a function of theta
Vcyl = pi * (B/2)^2 * d + V_c; % [m^3] Free cylinder volume as a function of theta

P = ones(1, 231) * P_atm; % Constant pressure throughout intake stroke