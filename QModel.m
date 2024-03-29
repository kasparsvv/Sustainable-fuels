function [Q] = QModel(Ca(i),CaS,CaD,mfuel,Q LHV g);
%Qmodel(Ca):: computes heat release by combustion
a = 5;  
%   Input: Ca, crank angle
global Runiv
if (isempty(Runiv))
    fprintf('[Qmodel] Assign global Runiv\n');
    return
end
xb(i) = 1-exp(-a*((Ca(i)-CaS)/CaD)^n);                                                  %Amount of fuel converted, Formula from project handbook page 13
dQcomb_di(i) = Q_LHV_E0 * mfuel * n * a * ((1-xb(i))/CaD) * ((Ca(i)-CaS)/CaD)^(n-1);    %Heat release, Formula from project handbook page 14
Q=dQcomb_di(i);

end
