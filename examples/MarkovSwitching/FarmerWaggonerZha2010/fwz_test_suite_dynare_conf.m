%% housekeeping
close all
clear
clc
%% add the necessary paths
setpaths

%% create the generic model object
generic=rise('fwz10_nasty');
%% 3 parameterizations
% the first parameterization has a unique solution
pp_1={'name','value','regime'
    'delta',    0,      1
    'delta',    0,      2
    'betta',    1,      1
    'betta',    1,      2
    'rho',     .9,      1
    'rho',     .9,      2
    'phi',     .5,      1
    'phi',     .8,      2
    'a_tp_1_2',1-.8,    1
    'a_tp_2_1',1-.9,    1};   

% the second parameterization has 2 solutions
pp_2={'name','value','regime'
    'delta',  -.7,       1
    'delta',   .4,       2
    'betta',    1,       1
    'betta',    1,       2
    'rho',      0,       1
    'rho',      0,       2
    'phi',     .5,       1
    'phi',     .8,       2
    'a_tp_1_2', 1-1,     1
    'a_tp_2_1',1-.64,    1};   

% the third parameterization has more than 2 solutions
pp_3={'name','value','regime'
    'delta',   -.7,      1
    'delta',   -.2,      2
    'betta',     1,      1
    'betta',     1,      2
    'rho',       0,      1
    'rho',       0,      2
    'phi',      .2,      1
    'phi',      .4,      2
    'a_tp_1_2',1-.9,     1
    'a_tp_2_1',1-.8,     1};   

%% vector of models with the different parameterizations
solvers={'newton_kronecker','newton_kronecker_iteration','functional_iteration'};
number_of_solvers=numel(solvers);

for ii=1:3
    if ii==1
        eval(['models_with_solver_',int2str(ii),'=rise.empty(0);'])
    end
    pp_ii=eval(['pp_',int2str(ii)]);
    for solver=1:number_of_solvers
        eval(['models_with_solver_',int2str(solver),'(ii,1)=generic.set_parameters(pp_ii);'])
        eval(['models_with_solver_',int2str(solver),...
            '(ii,1)=models_with_solver_',int2str(solver),'(ii,1).set_options(''solver'',solvers{solver});'])
    end
end

%% solve the models for each parameterization conditional on starting at 0

for ii=1:3
    keyboard
    clc
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