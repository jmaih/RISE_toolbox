%% Housekeeping
clear
clc
close all

%% RISE the model

m=rise('nk','max_deriv_order',2);

%% Solve the model using the Maih perturbation

mm=solve(m);

%% solve the model using the FRWZ's naïve perturbation

% mu is a switching parameter and enters the steady state but it is ignored
% and all the switching parameters are perturbed
mfrwz1=solve(m,'solve_perturbation_type','frwz');

%% solve the model using the FRWZ perturbation

% mu is a switching parameter and enters the steady state and is the only
% perturbed parameter
mfrwz2=solve(m,'solve_perturbation_type',{'frwz',{'mu'}});

%% FRWZ perturbation + Groebner basis
% When it does not solve, well, it does not find even a single solution

mfrwz3=solve(m,'solve_perturbation_type',{'frwz',{'mu'}},...
    'solver','dsge_groebner','solve_check_stability',false);

%% solve the model using the Maih-Waggoner perturbation

mmw=solve(m,'solve_perturbation_type','mw','solver','dsge_schur');

%% Find all solutions of the Maih-Waggoner perturbation using dsge_schur

refine=true;checkStab=false;allSols=true;msvOnly=true;xplosRoots=false;debug=false;

mmw2=solve(m,'solve_perturbation_type','mw',...
    'solver',{'dsge_schur',refine,checkStab,allSols,msvOnly,xplosRoots,debug});

%% Find all solutions of the Maih-Waggoner perturbation using dsge_udc
clc
refine=true;checkStab=false;allSols=true;msvOnly=true;xplosRoots=false;debug=false;

mmw3=solve(m,'solve_perturbation_type','mw',...
    'solver',{'dsge_udc',refine,checkStab,allSols,msvOnly,xplosRoots,debug});

%% print the solution

print_solution(mm)

print_solution(mfrwz1)

print_solution(mfrwz2)

print_solution(mmw)

print_solution(mmw2)
