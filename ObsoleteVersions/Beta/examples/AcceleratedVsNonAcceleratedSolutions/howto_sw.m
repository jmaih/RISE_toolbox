%% housekeeping
close all
clear
clc
% this files needs updating and would not work with the current version of RISE
%% add the necessary paths
setpaths
%% construct 2 models objects:
% the first model is with the accelerated solution(default)
% sometimes, many iterations are required for finding the loose commitment
% solution. So we jack up the maximum number of iterations to 5000
sw_com=model('usmodel_ms_lc','MaxIter',5000);
% we also set the model to full commitment
sw_com=sw_com.set_parameters(struct('name','ProbabilityOfCommitment','value',1));

% the second model is without acceleration
sw_com(2)=sw_com(1).set_options('accelerate_solver',false);
%% list of variables of interest

varlist=char('mc','zcap','rk','k','pk','c','inve','y','lab','pinf','w','r');

%% now we solve the models under full commitment
t0=zeros(1,2);
for ii=1:2
    tic, sw_com(ii)=sw_com(ii).solve; t0(ii)=toc;
    sw_com(ii).print_solution(varlist)
end
disp('comparison of accelerated(1) vs non-accelerated(2) solutions under full commitment')
disp(t0/min(t0))

disp('N.B: the commitment solutions will be different in terms of the')
disp('multipliers. This is probably due to the fact that the matrix we apply')
disp('the qr factorization on is not full. But, beyond and above everything,')
disp('what is interesting is that the solutions are 100 percent equivalent')
%% we repeat the exercise above but this time under discretion
sw_discr=sw_com;
sw_discr(1)=sw_discr(1).set_parameters(struct('ProbabilityOfCommitment',0));

sw_discr(2)=sw_discr(2).set_parameters(struct('ProbabilityOfCommitment',0));

t0=zeros(1,2);
for ii=1:2
    tic, sw_discr(ii)=sw_discr(ii).solve; t0(ii)=toc;
    sw_discr(ii).print_solution(varlist)
end
disp('comparison of accelerated(1) vs non-accelerated(2) solutions under discretion')
disp(t0/min(t0))
