%% Housekeeping
clear all
close all
clc
%% List of models and solvers
Models={'frwz_nk','frwz_nk@cal_2'};
solvers={'newton_system','newton_kronecker',...
    'newton_kronecker_iteration','functional_iteration'};
	
%% Alternative calibrations
% psi is controled by markov chain "a" and assumes a value of 2 in the second state
cal_2=struct();
cal_2.psi_a_2=0.7;

%% read the models and their calibrations
generic=cell(1,numel(Models));
for imod=1:numel(Models)
    mm=Models{imod};
    alphakrol=strfind(mm,'@');
    if isempty(alphakrol)
        generic{imod}=rise(mm);
    else
        generic{imod}=rise(mm(1:alphakrol-1));
        generic{imod}=set(generic{imod},'parameters',eval(mm(alphakrol+1:end))); 
    end
end
%% set the initialization scheme
backward=generic;
origin=generic;
for imod=1:numel(Models)
    backward{imod}=set_options(backward{imod},'solve_initialization','backward');
    origin{imod}=set_options(origin{imod},'solve_initialization','backward');
end
%% solve using all solvers
Solutions_backward=cell(numel(Models),numel(solvers));
Solutions_origin=cell(numel(Models),numel(solvers));
SolutionTimes_backward=zeros(numel(Models),numel(solvers));
SolutionTimes_origin=zeros(numel(Models),numel(solvers));
Nsim=500;
for isol=1:numel(solvers)
    for isim=1:Nsim
        tic,Solutions{1,isol}=solve(bk,[],'solver',solvers{isol});SolutionTimes(1,isol)=SolutionTimes(1,isol)+toc;
        tic,Solutions{2,isol}=solve(zer,[],'solver',solvers{isol});SolutionTimes(2,isol)=SolutionTimes(2,isol)+toc;
    end
end

SolutionTimes=SolutionTimes/Nsim; 
% conclusions:
% 1- the backward initialization is better for all the Newton algorithms, it
% does not matter for functional iteration.
% 2- the kronecker newton algorithms are the fastest, most certainly
% because the model is small. Functional iteration comes second in terms of
% computational speed.
% 3- the gains of using the Tadonki scheme are not visible here because:
% a) the model is small and the overhead from calling the extra routines
% for computing the TFQMR is non-negligeable
% b) the programming of the TFQMR and the kronecker times vector could
% probably be improved.

%     0.0112    0.0112    0.0207    0.0169
%     0.0114    0.0114    0.0226    0.0169
%% changing the baseline parameterization: 
% this new parameterization has two stable solutions and of course the
% question is, which one will be picked...
bk2=set(bk,'parameters',cal_2); 
zer2=set(zer,'parameters',cal_2);
for isol=1:numel(solvers)
    disp([' =============== ',solvers{isol},' =============== '])
    sbk=solve(bk2,'solver',solvers{isol});
    szer=solve(zer2,'solver',solvers{isol});
    sbk.print_solution
    szer.print_solution
end

% the conclusion is that 
% 1- all algorithms still find the same solution. In a sense then the zero
% initialization seems consistent with the backward initialization as they
% produce the same solution.  
% 2- But we know that there is another stable solution. And the question is
% why this particular solution has been chosen and is consistent with both
% zero and backward initializations? Is it the most stable of the
% solutions?
%% find all the solutions for the second parameterization.
% In this case, initialization does not matter any more
AllSolutions=cell(1,numel(solvers));
for isol=1:numel(solvers)
    disp([' =============== ',solvers{isol},' =============== '])
    newmod=set_options(frwz_nk,'solver',solvers{isol});
    newmod=set(newmod,'parameters',cal_2);
    AllSolutions{isol}=newmod.solve_alternatives;
end
% With random initialization, the Newton methods find both solutions.
% Functional iteration only finds one solution. It must be the case that
% the alternative solution violates some of the conditions of the theorem.
%% check whether there exist initial guesses that solve
% Now, given that we know the solution that is not found, the question is
% whether we can locate a starting value for which functional iteration
% would solve.

[criterion,aplus,aminus,tau,delta_zero,delta_bk]=...
    functional_iteration_convergence_conditions(...
    AllSolutions{1}(1));
disp(criterion)


