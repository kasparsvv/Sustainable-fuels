clear all;close all;clc
%% add general to matlab path 
addpath('General');
%%
DataDir='../Data/Training Set';  % The directory with the files. (Make sure you have backups!)
% DataDir='../Data/Gasoline';  % The directory with the files. (Make sure you have backups!)
%%
Files=dir(fullfile(DataDir,'*.txt'));nFiles=length(Files);                  % dir gives a directory listing, only *.txt files in this case
h=waitbar(0,'Converting');
for i=1:nFiles
    waitbar(i/nFiles,h);
    fname       = Files(i).name;                                            % Take a name from the list
    curfilename = fullfile(DataDir,fname);                                  % Create the full name
    comma2point_overwrite( curfilename);    % Convert comma to dot decimal separator, this is extremly fast!!
end
delete(h)
