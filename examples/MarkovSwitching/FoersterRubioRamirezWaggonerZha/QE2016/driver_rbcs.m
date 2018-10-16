%% Housekeeping
clear
clc
close all

%% RISE the model

m=rise('rbcs','max_deriv_order',2);

%% solve the model using the FRWZ perturbation approach

% mu is a switching parameter and enters the steady state
% steady-state partitioning
msp=solve(m,'solve_perturbation_type',{'frwz',{'mu'}});

% naive pertubation
mn=solve(m,'solve_perturbation_type','frwz');

%% print the solution

print_solution(msp)

print_solution(mn)
% N.B: The Naive perturbation solution in the paper is wrong

%% Solve for the second-order perturbation

msp2=solve(msp,'solve_order',2);

mn2=solve(mn,'solve_order',2);

%% print the solution

print_solution(msp2)

print_solution(mn2)