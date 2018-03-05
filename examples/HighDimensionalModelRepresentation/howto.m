%% housekeeping
clearvars
close all
clc

%% choose some options
% scramble the pseudo-random numbers or not
scramble_flag=true;
% optimal polynomial or theoretical ones
optimal_poly=true;
% debug or not
debug_flag=true;
%% read information about the function of interest
[objective,bounds]=ishigami();
% [objective,bounds]=satelli_sobol95();
lb=bounds(:,1);
ub=bounds(:,2);
npar=numel(lb);
param_names=strcat('pp_',num2str((1:npar)'));
%% choose the polynomial order
pol_order=4;
%% choose the expansion_order
expansion_order=3;
%% choose the number of samples
N=2^9;%2048;
%% option 1: generate samples and create the hdmr object

% theta=quasi_monte_carlo.sobol(npar,N,lb,ub,scramble_flag); % number of sub-intervals
theta=quasi_monte_carlo.sobol(lb,ub,N,scramble_flag); % number of sub-intervals
f=objective(theta);
objective_={f,theta};
%% option 2: let the hdmr object generate the samples for you
% objective_={objective,bounds,N};
%% construct the object
obj=hdmr(objective_,param_names,bounds,expansion_order,pol_order,optimal_poly);

%% estimate the object
profile off
profile on
obj=estimate(obj,debug_flag);
profile off
profile viewer
%% plot the fit insample
plot_fit(obj,'insample');
%% plot the fit out of sample
plot_fit(obj,'outofsample');
%% plot the individual effects
for ii=1:size(theta,1)
    figure();
    first_order_effect(obj,'insample',ii);
end
