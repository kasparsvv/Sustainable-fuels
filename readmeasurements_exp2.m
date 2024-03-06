clear all;close all;clc
%% add general to matlab path 
addpath('General');
%%
DataDir='';ColumnOrder={'time','Sensor','Encoder'}; %here the directory should be to a folder with just one measurement for a fuel mixture
% at half no and full load
% DataDir='../Data/Gasoline';cOrder={'time'','Encoder','Sensor};
col = lines(3);
%% Measurements
    S=0.055;
    B=0.067;
    l=0.084;
    rc=8.5;
    r=S/2;
%% Loading all measurments in DataDir
figure(1)
Files=dir(fullfile(DataDir,'*.txt'));nFiles=length(Files);                  % dir gives a directory listing, only *.txt files in this case
for k=1:nFiles
    fname       = Files(k).name;                                            % Take a name from the list
    curfilename = fullfile(DataDir,fname);                                  % Create the full name
    Data        = ImportData4GB10(curfilename,ColumnOrder);                             % Read the data. Type help ImportData4GB10
    % store it for each case. Yes a struct that contains a struct and other
    % things. Why not?
    Case(k).Data     = Data;
    Case(k).filename = fname;
    Case(k).DataDir  = DataDir;
    % preamble, put data in easy to use arrays
    t      = Data.t; %time array
    p      = Data.pulse; %pulse array
    V      = Data.Volt; %sensor voltage array             
    RevEnd = Data.RevEnds; %indices to revolution end point
    revnr = Data.NRevs; %number of complete revolutions
%%
    pos1 = RevEnd(1);                                                   
    pos2 = RevEnd(3); %a complete cycle includes 2 revolutions
    %V_cycle = V(pos1:pos2);                                           
    %p_max = max(V_cycle);  %find max pressure value in the cycle 

    %position_max_pressurei = find(V_cycle==p_max);        
    %position_max_pressure = position_max_pressurei(1);
    
    %Doubletooth = deg2rad(16) - (position_max_pressure-1)*(2*pi/((pos2 - pos1+1)/2));  
    Doubletooth=deg2rad(-223); %put in double tooth angle
    pos1 = RevEnd(1);                                                   
    pos2 = RevEnd(3); %a complete cycle includes 2 revolutions

    Ca(1) = Doubletooth; %first crank angle is at double tooth
    Ca(RevEnd(1)) = Doubletooth; %crank angle at the end of a revolution is also at the double tooth
    dCa = 2*pi/((pos2-pos1+1)/2); %how much the crank angle varies per index 
    Ca(RevEnd(1)-1) = Doubletooth - dCa;
    Volume(RevEnd(1)) = Vcyl(Ca(1),S,B,l,rc); %volume at the end of a revolution is found using the Vcyl function
    Pressure(RevEnd(1)) = ((V(RevEnd(1)))/5-0.115)/0.00385;
   
    figure(k)  
    for j=3:2:(length(RevEnd)-1)                                                
        for i=(RevEnd(j-2)+1):RevEnd(length(RevEnd))                         
            if i~=RevEnd(j)                                                  
            Ca(i) = Ca(i-1) + dCa;
            Volume(i) = Vcyl(Ca(i),S,B,l,rc);
            Pressure(i) = ((V(i))/5-0.115)/0.00385;
            else
            Volume = Volume((RevEnd(j-2)):(RevEnd(j)-1));
            Pressure = Pressure((RevEnd(j-2)):(RevEnd(j)-1));
            plot(Volume,Pressure)
            xlabel('Volume [m^3]')
            ylabel('Pressure [bar]')
            legend('Real cycle')
            if k==1
                title('No load')
            else 
                if k==2
                    title('Full load')
            else title('Half load')
                end    
            end
            hold on
            Ca(i) = Ca(1);
            Volume = [];
            Pressure = [];
            Volume(i) = Vcyl(Ca(i),S,B,l,rc);
            Pressure(i) = ((V(i))/5-0.115)/0.00385;
            break
            end
        end
    end
end
%%
function V_cyl = Vcyl(Ca,S,B,l,rc)

r=S/2;
V_d=pi*(B/2)^2*S;
V_c=-V_d/(1-rc);

x = r*cos(Ca)+sqrt(l^2-r^2*sin(Ca)^2);
d=l+r-x;
V_cyl=pi*(B/2)^2*d+V_c;
end