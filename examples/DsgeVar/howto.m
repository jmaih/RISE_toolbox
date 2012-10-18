%% housekeeping
close all
clear
clc
%% add the necessary paths
setpaths
%% estimation of the model with

datarabanal_hybrid
data=rise_time_series(1,[pie r rw y],char('pie','r','rw','y'));
data=data.window(50,50+90-1);

mx=rise('bvar_forward_ms','data',data);

mx=mx.estimate(1); % 261.162 

%% simulate the posterior distribution

mx=mx.posterior_simulator(true);

% a sub-folder is created with the name, simulation_folder

%% compute IRFs: 
% 2 types of irfs are computed: the dsge irfs are stored in the usual way
% while the bvar_dsge irfs are stored under mx.obj.dsge_var

mx=mx.set_options('irf_draws',20);

% irfs drawing around the mode
mx_mode=mx.irf(40,1,[],'mode');

% irfs drawing from the posterior distribution
mx_post=mx.irf(40,1,[],'posterior');

