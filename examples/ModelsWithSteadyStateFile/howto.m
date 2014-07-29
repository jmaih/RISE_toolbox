%% housekeeping
close all
clear
clc
%% add the necessary paths
rise_startup()
%% read the model and add the name of the steady state file
mwssf=rise('Canonical_Const','steady_state_file','steady_state_4_Canonical_Const');

%% solve the model
mwssf=mwssf.solve;

%% SEE ALSO: NonlinearModels