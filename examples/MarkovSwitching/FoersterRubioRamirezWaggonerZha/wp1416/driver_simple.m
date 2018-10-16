%% Housekeeping
clear
clc
close all

%% RISE the model

m=rise('simple');

%% solve the model using the Naïve FRWZ perturbation approach
% phi and sigma do not affect the steady state
m1=solve(m,'solve_perturbation_type',{'frwz',{'phi','sigma'}});

%% solve the model using the FRWZ perturbation approach

m2=solve(m,'solve_perturbation_type','frwz');

%% print the solution

print_solution(m1)

print_solution(m2)