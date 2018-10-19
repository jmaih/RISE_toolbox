%% housekeeping
close all
clear
clc

%% create the generic model object
m0=rise('fwz10_nasty');
%% 3 parameterizations
% the first parameterization has a unique solution
pp=struct();
pp.delta_a_1=0; 
pp.delta_a_2=0; 
pp.betta_a_1=1; 
pp.betta_a_2=1; 
pp.rho_a_1=.9;
pp.rho_a_2=.9;
pp.phi_a_1=.5; % controlled by markov chain 'a' and assumes .5 in state 1 and .8 in state 2
pp.phi_a_2=.8;
pp.a_tp_1_2=1-.8; % for all regimes
pp.a_tp_2_1=1-.9; % for all regimes   

% the second parameterization has 2 solutions
pp(2)=pp;
pp(2).delta_a_1=-.7;
pp(2).delta_a_2=.4; 
pp(2).rho_a_1=0;
pp(2).rho_a_2=0;	
pp(2).a_tp_1_2=1-1;
pp(2).a_tp_2_1=1-.64;   

% the third parameterization has more than 2 solutions
pp(3)=pp(2);
pp(3).delta_a_2=-.2;
pp(3).phi_a_1=.2;
pp(3).phi_a_2=.4;
pp(3).a_tp_1_2=1-.9;
pp(3).a_tp_2_1=1-.8;  

np=numel(pp);

%% solve models with the different parameterizations conditional on starting at the backward solution
solvers={'mnk','fwz','mn','mfi'};

number_of_solvers=numel(solvers);

m=struct();

for ii=1:np
        
    for solver=1:number_of_solvers
        
        m(ii).(solvers{solver})=solve(m0,'parameters',pp(ii),...
            'solver',solvers{solver});
        
        disp(['Parameterization ::',int2str(ii),', solver ::',solvers{solver}])
        
        print_solution(m(ii).(solvers{solver}))
        
    end
    
end

%% for each parameterization, find all possible solutions
clc

M=struct();

for ii=1:numel(m)
    
    for solver=1:number_of_solvers
        
        M(ii).(solvers{solver})=solve_alternatives(m(ii).(solvers{solver}));
        
        disp(['Parameterization ::',int2str(ii),', solver ::',solvers{solver},...
            ' # solutions ',int2str(numel( M(ii).(solvers{solver})))])

    end
    
end

%% latest solvers
clc
refine=false;checkStab=false;allSols=true;msvOnly=true;
xplosRoots=false;debug=false;
solver_noref={'dsge_udc',refine,checkStab,allSols,msvOnly,xplosRoots,debug};
solver_ref={'dsge_udc',true,checkStab,allSols,msvOnly,xplosRoots,debug};

M2=cell(1,np);

for ii=1:np
    
    m0i=set(m0,'parameters',pp(ii),'solve_check_stability',false);
        
    m0i_mw=set(m0i,'solve_perturbation_type','mw');
    
    mi=rise.empty(0,1);
    % Groebner: original model
    %--------------------------
    mi(1)=solve(m0i,'solver','dsge_groebner');
    % Groebner: mw-perturbed model
    %------------------------------
    mi(2)=solve(m0i_mw,'solver','dsge_groebner');
    % dsge_udc: mw-perturbed model & no refinement
    %----------------------------------------------
    mi(3)=solve(m0i_mw,'solver',solver_noref);
    % dsge_udc: mw-perturbed model & with refinement
    %------------------------------------------------
    mi(4)=solve(m0i_mw,'solver',solver_ref);
    
    M2{ii}=mi;
end

