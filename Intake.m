%% Defining constants

S = 0.055; % [m] Stroke
r  = 1/2 * S; % [m] Lenght of crankshaft
l = 0.0842; % [m] Lenght of the connecting rod
B = 0.0677; % [m] Bore dimension
TDC = 0.003; % [m] Top dead center
BDC = S + TDC; % [m] Bottom dead center
V_d = (pi/4) * B^2 * S; % [m^3] Displacement volume
V_c = V_d/(rc - 1); % [m^3] Clearance volume
P_atm  = 1; % [Bar] Atmospheric pressure (assumed for now)

%% Defining the kinematic equations of the engine as functions of theta

theta_deg = linspace(0, 180, 181); % Intake stroke in degrees. Lasts from -15 to 215
theta_rad = deg2rad(theta_deg); % Intake stroke in radians for ease of use

x = r * cos(theta_rad) + sqrt(l^2 - r^2 * sin(theta_rad).^2); % [m] x postion of piston as a function of theta
d = l + r - x; % [m] Distance of pistion from TDC as a function of theta
Vcyl = pi * (B/2)^2 * d + V_c; % [m^3] Free cylinder volume as a function of theta

P = ones(1, 181) * P_atm; % Constant pressure throughout intake stroke

%% Plotting different figures

% PV
figure;
plot (Vcyl, P) ;
xlabel ( ' Volume ( m ^3) ') ;
ylabel ( ' Pressure ( Bar ) ') ;
title ( 'Intake stroke') ;
grid on ;

% Theta V
figure;
plot (theta_deg, Vcyl) ;
xlabel ( '\theta (degrees)') ;
ylabel ( 'Specific volume [m^3]') ;
title ( 'Intake stroke') ;
grid on ;

% Theta x
figure;
plot (theta_deg, x) ;
xlabel ( '\theta (degrees)') ;
ylabel ( ' x [m] ') ;
title ( 'Intake stroke') ;
grid on ;