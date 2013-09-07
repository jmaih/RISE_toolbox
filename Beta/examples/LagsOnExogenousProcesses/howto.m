%% housekeeping
close all
clear
clc
%% add the necessary paths
setpaths

%% Nothing special here. Look at the model file.

loe=rise('Canonical_Const');

loe=loe.solve;

