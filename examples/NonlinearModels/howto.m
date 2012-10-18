%% housekeeping
close all
clear
clc
%% add the necessary paths
disp(upper('this file may need updating...'))

setpaths
%% we pass a function that solves the steady state analytically

mwssf_a=rise('fs2000','steady_state_file','fs2000_steadystate');

mwssf_test=mwssf_a.solve;

profile off
profile on
mwssf_a=mwssf_a.solve;
profile off
profile viewer

%% we pass a function that gives initial values of the steady state. 
% this is the counterpart to dynare's initval

mwssf_b=rise('fs2000','steady_state_file','fs2000_steadystate_initval');

profile off
profile on
mwssf_b=mwssf_b.solve;
profile off
profile viewer

%% compare the solutions

max(max(abs(mwssf_a.T-mwssf_b.T)))

%% compare the steady states

[vertcat(mwssf_a.varendo.det_steady_state)-vertcat(mwssf_b.varendo.det_steady_state)]
