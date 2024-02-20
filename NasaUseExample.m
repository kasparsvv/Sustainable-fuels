%% Test data base
clear all;close all;
%% This global is needed by the Cp, Cv etc functions
addpath('General/Nasa'); % Add directory of Nasa routines to Matlab-path
global Runiv 
Runiv = 8.314;
%% 
DBdir = 'General\Nasa';
DBname = 'NasaThermalDatabase';
load(fullfile(DBdir,DBname));
Sp
El
%% I want to find all alcohols in the database
iCHO=myfind({El.Name},{'C','H','O'}); % Which position corresponds to C H and O elements in the Sp.Elcomp vector
iH=iCHO(2);iC=iCHO(1);iO=iCHO(3);
iAlc=[];ii=0;
for i=1:length(Sp)
    if ( Sp(i).Elcomp(iH)==(2*Sp(i).Elcomp(iC)+2) && Sp(i).Elcomp(iO)==1) % Typical comp of an alcohol 2*nC+2 H-atoms, 1 O-atom
        fprintf('%20s [%3i] \n',Sp(i).Name,i);
        ii = ii+1;
        iAlc(ii)=i;
    end;
end
%% Plot all Cp values for given T range for the found species
T=300:10:3000;nAlc=ii;
for j=1:nAlc
    Cp(j,:)=CpNasa(T,Sp(iAlc(j)));
    Cv(j,:)=CvNasa(T,Sp(iAlc(j)));
end
figure(1)
subplot(1,2,1)
plot(T,Cp,'-',T,Cv,'--','LineWidth',2);
legend(Sp(iAlc).Name);
xlabel('T [K]');ylabel('c_p [J/kg/K]');
subplot(1,2,2)
plot(T,Cp./Cv,'-','LineWidth',2);
legend(Sp(iAlc).Name);
xlabel('T [K]');ylabel('\gamma [-]');
%test