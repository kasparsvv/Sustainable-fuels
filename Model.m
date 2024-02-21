relativepath_to_generalfolder='General'; % relative reference to General folder (assumes the folder is in you working folder)
addpath(relativepath_to_generalfolder); 
% Load Nasadatabase
TdataBase=fullfile('General','Nasa','NasaThermalDatabase');
load(TdataBase);


%Constants
P0=1.01235e5 ;      % Ambient Pressure (pa)
T0=293;             % Ambient Temparature (K)
S=1;                % Stroke (m)
B=1;                % Bore (m)
l=1;                % Length of the connecting rod (m)
rc=1;               % Compression ratio (-)
Runiv=8.314472;     % Universal gas constant (J / mol·K)
CaS=1;              % Crank angle at start of combustion
CaD=1;              % Combustion duration (in terms of crank angle)
n=3;                % Wiebe form factor, Project handbook says 3 is often used
a=5;                % Wiebe efficiency factor, Project handbook says 5 is often used


% Select species for the case at hand
cFuel = 'Gasoline';
iSp = myfind({Sp.Name},{cFuel,'O2','CO2','H2O','N2'});                      % Find indexes of these species
SpS=Sp(iSp);                                                                % Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSp = length(SpS);
Mi = [SpS.Mass];
% Fuel composition
Xfuel = [0 0 0 0 0];                                                        % Order is important. Note that these are molefractions
Mfuel = Xfuel*Mi';                                                          % Row times Column = inner product 
Yfuel = Xfuel.*Mi/Mfuel;                                                    % Vector. times vector is Matlab's way of making an elementwise multiplication

Rg = Runiv/Mfuel;   %Specific gas constant


% Initialisation
p(1)=P0;T(1)=T0;
pad(1)=p(1);
Ca(1)=0.0;
V(1)=Vcyl(Ca(1),S,B,l,rc); % Vcyl is a function that computes cylinder volume as function of crank-angle for given B,S,l and rc
m(1)=p(1)*V(1)/Rg/T(1);


% Loop over crank-angle, with 'for' construction
NCa=360;            % Number of crank-angles
dCa=0.5;            % Stepsize
NSteps=NCa/dCa;

%For loop to calculate e.g. volume at given crank-angle
for i=2:NSteps,                         % Calculate values for 1 cycle
    Ca(i)=Ca(i-1)+dCa;
    V(i)=Vcyl(Ca(i),S,B,l,rc);          % New volume for current crank-angle
    m(i)=m(i-1);                        % Mass is constant, valves are closed
    dV=V(i)-V(i-1);                     % Volume change
    dQcom = QModel(Ca(i));              % Heat Release by combustion
    dT=(dQcom-p(i-1)*dV)/Cv/m(i-1);     % 1st Law dU=dQ-pdV (closed system)
                                        % adiabatic closed system with constant
                                        % gas composition and constant Cv
    T(i)=T(i-1)+dT;
    p(i)=m(i)*Rg*T(i)/V(i);             % Gaslaw


    for n=1:NSp
        Cvi(n) =HNasa(T(i),SpS(n));           % Get Cv from Nasa-table
    end
    %Cvi=YourNasaThermofunction(T,Sp); % Nasa polynomials (see JetEngine)


    Cv=sum(Yi.*Cvi); % heat cap at current T
end;
