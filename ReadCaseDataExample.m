clear all;close all;clc
%% add general to matlab path 
addpath('General');
%%
DataDir='../Data/Training Set';ColumnOrder={'time','Sensor','Encoder'};
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
    V      = Data.Volt;
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