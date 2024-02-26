%% Defining constants

S = 0.055; % [m] Stroke
r  = 1/2 * S; % [m] Lenght of crankshaft
l = 0.0842; % [m] Lenght of the connecting rod
B = 0.0677; % [m] Bore dimension
TDC = 0.003; % [m] Top dead center
BDC = S + TDC; % [m] Bottom dead center
V_c = pi/4 * B^2 * TDC; % [m63] Compression volume (will change)
P_atm  = 1; % [Bar] Atmospheric pressure (assumed for now)

%% Defining the kinematic equations of the engine as functions of theta

theta_deg = linspace(0, 180, 181); % Intake stroke in degrees. Lasts from -15 to 215
theta_rad = deg2rad(theta_deg); % Intake stroke in radians for ease of use

x = r * cos(theta_rad) + sqrt(l^2 - r^2 * sin(theta_rad).^2); % [m] x postion of piston as a function of theta
d = l + r - x; % [m] Distance of pistion from TDC as a function of theta
Vcyl = pi * (B/2)^2 * d + V_c; % [m^3] Free cylinder volume as a function of theta

P = ones(1, 181) * P_atm; % Constant pressure throughout intake stroke

%% Plotting different figures

% % PV
% figure;
% plot (Vcyl, P) ;
% xlabel ( ' Volume ( m ^3) ') ;
% ylabel ( ' Pressure ( Bar ) ') ;
% title ( 'Intake stroke') ;
% grid on ;
% 
% % Theta V
% figure;
% plot (theta_deg, Vcyl) ;
% xlabel ( '\theta (degrees)') ;
% ylabel ( 'Specific volume [m^3]') ;
% title ( 'Intake stroke') ;
% grid on ;
% 
% % Theta x
% figure;
% plot (theta_deg, x) ;
% xlabel ( '\theta (degrees)') ;
% ylabel ( ' x [m] ') ;
% title ( 'Intake stroke') ;
% grid on ;


% Compression stroke, ignition and Power Stroke

C1 = p(1) * V(1) ^ gamma; % Poisson relations
C2 = T(1) * V(1) ^(gamma - 1); % Poisson relations

% Compression
if Ca(i) >= -180 && Ca(i) < 0
    p(i) = C1 / V(i) ^(gamma); % Poisson relations
    T(i) = C2 / V(i) ^(gamma - 1); % Poisson relations
end

% Ignition
if Ca(i) == 0
    dQcom(i) = 730; % Heat Release by combustion
    dT(i) = (dQcom(i) - p(i - 1) * dV) / Cv / m(i - 1); % 1st Law dU=dQ -pdV (closed system)
    T(i) = T(i - 1) + dT(i);
    p(i) = m(i) * Rg * T(i) / V(i); % Gaslaw
    C3 = p(i) * V(i) ^ gamma; % Poisson relations
    C4 = T(i) * V(i) ^(gamma - 1); % Poisson relations
end

% Power stroke
if Ca(i) > 0 && Ca(i) < 180
    p(i) = C3 / V(i) ^(gamma); % Poisson relations
    T(i) = C4 / V(i) ^(gamma - 1); % Poisson relations
end

