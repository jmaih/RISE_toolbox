%% housekeeping
close all
clear
clc
%% add the necessary paths
rise_startup()
%% set the profiler
profile off
profile on
%% create the generic model object
generic=rise('fwz10_nasty');
%% 3 parameterizations
% the first parameterization has a unique solution
pp_1=struct();
pp_1.delta=0; % assumes a value of 0 in all states
pp_1.betta=1; % assumes a value of 1 in all states
pp_1.rho=.9;
pp_1.phi_a_1=.5; % controled by markov chain 'a' and assumes .5 in state 1 and .8 in state 2
pp_1.phi_a_2=.8;
pp_1.a_tp_1_2=1-.8; % for all regimes
pp_1.a_tp_2_1=1-.9; % for all regimes   

% the second parameterization has 2 solutions
pp_2=struct();
pp_2.delta_a_1=-.7;
pp_2.delta_a_2=.4;
pp_2.betta=1;
pp_2.rho=0;
pp_2.phi_a_1=.5;
pp_2.phi_a_2=.8;
pp_2.a_tp_1_2=1-1;
pp_2.a_tp_2_1=1-.64;   

% the third parameterization has more than 2 solutions
pp_3=struct();
pp_3.delta_a_1=-.7;
pp_3.delta_a_2=-.2;
pp_3.betta=1;
pp_3.rho=0;
pp_3.phi_a_1=.2;
pp_3.phi_a_2=.4;
pp_3.a_tp_1_2=1-.9;
pp_3.a_tp_2_1=1-.8;   

%% vector of models with the different parameterizations
solvers={'newton_kronecker','newton_system','newton_kronecker_iteration',...
    'functional_iteration'};
number_of_solvers=numel(solvers);

for ii=1:3
    if ii==1
        eval(['models_with_solver_',int2str(ii),'=rise.empty(0);'])
    end
    pp_ii=eval(['pp_',int2str(ii)]);
    for solver=1:number_of_solvers
        eval(['models_with_solver_',int2str(solver),'(ii,1)=generic.set(''parameters'',pp_ii);'])
        eval(['models_with_solver_',int2str(solver),...
            '(ii,1)=models_with_solver_',int2str(solver),'(ii,1).set_options(''solver'',solvers{solver});'])
    end
end

%% solve the models for each parameterization conditional on starting at 0

for ii=1:3
    for solver=1:number_of_solvers
        eval(['models_with_solver_',int2str(solver),'(ii,1)=models_with_solver_',int2str(solver),'(ii,1).solve;'])
        disp(['Parameterization ::',int2str(ii),', solver ::',solvers{solver}])
        eval(['models_with_solver_',int2str(solver),'(ii,1).print_solution'])
    end
end

%% for each parameterization, find all possible solutions
Final_Results=cell(4,number_of_solvers+1);
for ii=1:3
    Final_Results{ii+1,1}=['pp_',int2str(ii)];
    for solver=1:number_of_solvers
        if ii==1
            Final_Results{1,solver+1}=solvers{solver};
        end
        disp(['Parameterization ::',int2str(ii),', solver ::',solvers{solver}])
        eval(['bank_solver_param_',int2str(ii),'_',int2str(solver),'=models_with_solver_',int2str(solver),'(ii,1).solve_alternatives;'])
        Final_Results{ii+1,solver+1}=max(size(eval(['bank_solver_param_',int2str(ii),'_',int2str(solver)])));
    end
end
disp(Final_Results)
profile off 
profile viewer