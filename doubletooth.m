clear all;close all;clc
%% add general to matlab path 
addpath('General');
%%
DataDir='';ColumnOrder={'time','Sensor','Encoder'}; %here the directory should be to a folder with all the engine turning off measurements
% DataDir='../Data/Gasoline';cOrder={'time'','Encoder','Sensor};
col = lines(3);
%% Loading all measurments in DataDir
figure(1)
Files=dir(fullfile(DataDir,'*.txt'));nFiles=length(Files);                  % dir gives a directory listing, only *.txt files in this case
for i=1:nFiles
    fname       = Files(i).name;                                            % Take a name from the list
    curfilename = fullfile(DataDir,fname);                                  % Create the full name
    Data        = ImportData4GB10(curfilename,ColumnOrder);                             % Read the data. Type help ImportData4GB10
    % store it for each case. Yes a struct that contains a struct and other
    % things. Why not?
    Case(i).Data     = Data;
    Case(i).filename = fname;
    Case(i).DataDir  = DataDir;
    % preamble, put data in easy to use arrays
    t      = Data.t;
    p      = Data.pulse;
    V      = Data.Volt;                                                       % Maximum pressure at 14 degrees after TDC
    RevEnd = Data.RevEnds;
    % Actual Plotting use of subplot
    figure(1)
    subplot(nFiles,2,2*(i-1) + 1);                                          % help subplot if you don't understand
    plot(t,p,'b-');
    line(t(RevEnd),p(RevEnd),'LineStyle','none','Marker','s','MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor','y');
    xlabel('t [s]');ylabel('pulse [V]');
    legend('pulse','double tooth');
    title(fname);
    %
    subplot(nFiles,2,2*(i-1) + 2);
    plot(t,V,'r-');
    line(t(RevEnd),V(RevEnd),'LineStyle','none','Marker','s','MarkerSize',8,'MarkerEdgeColor','r','MarkerFaceColor','y');
    xlabel('t [s]');ylabel('p Sensor Signal [V]');
    legend('p signal','location of double tooth');
    title(fname);


end
%% All measurements are loaded
for i = 1:nFiles
        pos1 = find(t==coordinateofdoubletooth1); %put in the coordinate of a double tooth
        pos2 = find(t==coordinateofdoubletooth2);
        V_cut = V(pos1:pos2);
        p_cut = p(pos1:pos2);
        t_cut = t(pos1:pos2);

figure(3)
    subplot(2,1,1)
    plot(t_cut,p_cut,'b-');
    line(t(RevEnd),p(RevEnd),'LineStyle','none','Marker','s','MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor','y');
    xlabel('t [s]');ylabel('pulse [V]');
    legend('pulse','double tooth');
    title(fname);
    subplot(2,1,2)
    plot(t_cut,V_cut,'r-');
    line(t(RevEnd),V(RevEnd),'LineStyle','none','Marker','s','MarkerSize',8,'MarkerEdgeColor','r','MarkerFaceColor','y');
    xlabel('t [s]');ylabel('p Sensor Signal [V]');
    legend('p signal','location of double tooth');
    title(fname);

        Max_pressure = max(V_cut);
        position_max_pressurei = find(V_cut==Max_pressure);
        position_max_pressure = position_max_pressurei(1);
        Doubletooth = deg2rad(0) - (position_max_pressure-1)*(2*pi/((pos2-pos1)/2));
        Doubletooth = rad2deg(Doubletooth);
end
