%% Housekeeping
clc
clear
close all

%% Using the Grobner solver: constant-parameter case
m=rise('rbcs0');
v0=solve(m);
print_solution(v0)

% Finding all solutions
%-----------------------
v1=solve(m,'solver','dsge_groebner','solve_check_stability',false);
print_solution(v1)
print_solution(v0)

%% Using the Grobner solver: regime switching example
clc

m=rise('rbcs','solve_perturbation_type',{'frwz',{'mu'}});

% Finding one solution
%----------------------
v0=solve(m);

% Finding all solutions
%-----------------------
v1=solve(m,'solver','dsge_groebner','solve_check_stability',false);
print_solution(v0)
print_solution(v1)