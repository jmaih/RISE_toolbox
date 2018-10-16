%% Housekeeping
clear
clc
close all

%% RISE the model

m=rise('nk','max_deriv_order',2);

%% solve the model using the FRWZ perturbation approach

% mu is a switching parameter and enters the steady state

m1=solve(m,'solve_perturbation_type',{'frwz',{'mu'}});

%% solve the model using the Maih-Waggoner perturbation approach

m2=solve(m,'solve_perturbation_type','mw');

%% print the solution

print_solution(m1)

%% Solve for the second-order perturbation

m2=solve(m1,'solve_order',2);

%% print the solution

print_solution(m2)