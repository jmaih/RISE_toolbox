%% housekeeping
close all
clear
clc
%% add the necessary paths
rise_startup()
%% read the model 

test=rise('Canonical_osr');

%% estimate with fmincon
test=estimate(test);

