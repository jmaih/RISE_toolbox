%% Housekeeping
clear 
close all
clc
%%
tic
%% create data

run('tut01_data')

%% estimate models

run('tut02_estimation')

%% compare residuals

run('tut03_compare_residuals')

%% set identification

run('tut04_identification')

%% Compare structural shocks

run('tut04b_compare_structural_shocks')
%% simple choleski IRFs

run('tut05_compare_choleski_irfs')

%% more involved restrictions

run('tut06_compare_combos_irfs')

%% variance decomposition

run('tut07_variance_decomposition')

%% historical decomposition

run('tut08_historical_decomposition')

%% sample distribution of parameters

run('tut09_bootstrap')

%% distribution of variance decomposition

run('tut10_variance_decomposition_distribution')

%% distribution of historical decomposition

run('tut11_historical_decomposition_distribution')

%% distribution of IRFs

run('tut12_irf_distribution')

%% Bayesian estimation

run('tut13_bayestimation')

%% forecast

run('tut14_bayesforecast')

%% conditional forecast

run('tut15_conditionalforecast')

%% entropic tilting

%% clean up
toc

delete *.mat
