%% housekeeping
clear
close all
clc
%% RISE the model
m=rise('nkp');

%% calibrated parameters and priors

start_from_mode=false;

[p,priors]=create_parameters(start_from_mode);

m=set(m,'parameters',p);

%% data

[data]=create_data();

%% set options

clc
log_vars=[];%{'A','C','G','LAMBDA','PAI','Q','R','THETA','X','Y','Z'};
m=set(m,'data',data,...
    'solve_log_approx_vars',log_vars,...
    'estim_priors',priors);

%% estimate model without endogenous priors

mest=estimate(m);

%% estimate model with endogenous priors
clc
% mest2=estimate(m,'estim_endogenous_priors',@simple_endo_priors);

[mest2,filts]=estimate(m,'estim_endogenous_priors',@cov_endo_priors);

%% How well did we do?

cov_endo_priors(mest2,filts)

%% Impulse response functions
myirfs=irf([mest,mest2]);

tex=get(m,'tex');

close all
% mest=set(mest,'tex_name',{
%     'Y','Output'
%     'PAI','Inflation'
%     'X','Output gap'
%     'R','Interest rate'
%     'EPS_E','Cost push'
%     'EPS_R','monetary policy'
%     'EPS_Z','Technology'
%     'EPS_A','Preference'
%     });
figure('name','Impulse responses');
myvars={'Y','PAI','R','X'};
nvars=numel(myvars);
iter=0;
for ivar=1:nvars
    vname=myvars{ivar};
    for ishock=1:mest.exogenous.number(1)
        shock=mest.exogenous.name{ishock};
        iter=iter+1;
        subplot(4,4,iter)
        plot('0:20',100*myirfs.(shock).(vname),'linewidth',2)
        axis tight
        title([tex.(vname),' to ',tex.(shock)])
        if ivar==1 && ishock==1
            legend('plain','endogenous priors')
        end
    end
end
orient landscape
