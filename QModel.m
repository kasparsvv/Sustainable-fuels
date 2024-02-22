function [Q] = QModel(Ca)
%Qmodel(Ca):: computes heat release by combustion
% 
%   Input: Ca, crank angle
global Runiv
if (isempty(Runiv))
    fprintf('[Qmodel] Assign global Runiv\n');
    return
end
xb(i) = 1-exp(-a*((Ca(i)-CaS)/CaD)^n);                                              %Amount of fuel converted, Formula from project handbook page 13
dQcomb_di(i) = Qlhv * mfuel * n * a * ((1-xb(i))/CaD) * ((Ca(i)-CaS)/CaD)^(n-1);    %Heat release, Formula from project handbook page 14
Q=dQcomb_di(i);

end