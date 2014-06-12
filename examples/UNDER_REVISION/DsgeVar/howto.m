%% housekeeping
close all
clear
clc
%% add the necessary paths
rise_startup()
%% estimation of the model with

datarabanal_hybrid
data=ts(1,[pie r rw y],char('pie','r','rw','y'));
data=data.window(50,50+90-1);

mx=rise('bvar_forward_ms','data',data);

mx=mx.estimate(); % 261.162 

%% simulate the posterior distribution

mx=posterior_simulator(mx);

%% computation of Marginal data density
% the result will be stored back into the
% object

mx=log_marginal_data_density_mhm(mx);

mx=log_marginal_data_density_chib_jeliazkov(mx);

