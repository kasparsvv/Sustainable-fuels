clear all;
close all;
run("Model.m")


RPM = 3000; % rounds per minute
n_rpc = 2; %number of rounds per cycle

W_E0= trapz(V,p);
efficiency = W_E0/dQcom*100;
gamma_average = (gamma_comb_out+gamma_comp_in) / 2;
ottoefficiency =(1-(1/rc)^(gamma_average-1)) *100;
P_E0 = W_E0 * (RPM/60) * (1/n_rpc);
