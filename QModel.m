function [H,HU] = QModel(Ca)
%Qmodel(Ca):: computes heat release by combustion
% 
%   Input: Ca, crank angle
global Runiv
if (isempty(Runiv))
    fprintf('[Qmodel] Assign global Runiv\n');
    return
end
xb(i)= 1-exp(-a*((Ca(i)-CaS)/CaD)^n);          %Amount of fuel converted, Formula from project handbook page 13
%dQcomb/di = ;
end