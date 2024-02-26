%% Defining constants

S = 0.055; % [m] Stroke
r  = 1/2 * S; % [m] Length of crankshaft
l = 0.0842; % [m] Length of the connecting rod
B = 0.0677; % [m] Bore dimension
TDC = 0.003; % [m] Top dead center
BDC = S + TDC; % [m] Bottom dead center
V_d = (pi/4) * B^2 * S; % [m^3] Displacement volume
rc = 10; % Compression ratio (assumed for now)
V_c = V_d/(rc - 1); % [m^3] Clearance volume
P_atm  = 1; % [Bar] Atmospheric pressure (assumed for now)

%% Defining the kinematic equations of the engine as functions of theta

Ca = linspace(0, 180, 181); % Intake stroke in degrees.

for i = 1:length(Ca)
    x(i) = r * cosd(Ca(i)) + sqrt(l^2 - r^2 * sind(Ca(i))^2); % [m] x position of piston as a function of theta
    d(i) = l + r - x(i); % [m] Distance of piston from TDC as a function of theta
    Vcyl(i) = pi * (B/2)^2 * d(i) + V_c; % [m^3] Free cylinder volume as a function of theta
end

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