%% Housekeeping
clc
clear
close all

%% Using the Grobner solver: constant-parameter case
m0=rise('rbcs0');
m1=solve(m0);
print_solution(m1)

% Finding all solutions
%-----------------------
m2=solve(m0,'solver','dsge_groebner','solve_check_stability',false);
print_solution(m2)
print_solution(m1)

%% Using the Grobner solver: regime switching example
clc

m00=rise('rbcs','solve_perturbation_type',{'frwz',{'mu'}});

% Finding one solution
%----------------------
m01=solve(m00);

% Finding all solutions
%-----------------------
m02=solve(m00,'solver','dsge_groebner','solve_check_stability',false);
print_solution(m01)
print_solution(m02)


%% Using the dsge_udc solver: regime switching example
clc

m000=set(m00,'solve_perturbation_type','mw');

% Finding one solution
%----------------------
m001=solve(m000);

% Finding all solutions Groebner
%-------------------------------
m002=solve(m000,'solver','dsge_groebner','solve_check_stability',false);

% Finding all solutions UDC
%--------------------------
refine=false;checkStab=false;allSols=true;msvOnly=true;
xplosRoots=false;debug=false;
solver={'dsge_udc',refine,checkStab,allSols,msvOnly,xplosRoots,debug};
m003=solve(m000,'solver',solver,'solve_check_stability',false);

% Finding all solutions UDC with refinement
%------------------------------------------
refine=true;checkStab=false;allSols=true;msvOnly=true;
xplosRoots=false;debug=false;
solver={'dsge_udc',refine,checkStab,allSols,msvOnly,xplosRoots,debug};

m004=solve(m000,'solver',solver,'solve_check_stability',false);

