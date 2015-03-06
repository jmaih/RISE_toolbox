%% housekeeping
clc
%% RISE the model
m=rise('targets','steady_state_file','sstate_model');
%% get the parameters
[p,priors]=create_parameters('a',false);
%% push the parameters
m=set(m,'parameters',p);
%% create data
[data]=create_data();
%% estimate the model

ms=estimate(m,'data',data,'estim_priors',priors);
%%
[mtest,LogLik3,~,retcode]=filter(m,'data',data);
%% Impulse responses
myirfs=irf(ms);
%%
myvars={'Y','PAI','R','X'};
locs=locate_variables(myvars,m.endogenous.name);
myvtex=m.endogenous.tex_name(locs);
sstate=get(mtest,'sstate');
ssdev=false;

shocks={'EPS_A','EPS_E','EPS_Z','EPS_V','EPS_PAI'};%m.exogenous.name;
locs=locate_variables(shocks,m.exogenous.name);
myshtex=m.exogenous.tex_name(locs);
close all
figure('name','Impulse reponses');
iter=0;
ss=1;
for ishock=1:numel(shocks)
    shock=shocks{ishock};
    for ivar=1:numel(myvars)
        vname=myvars{ivar};
        if ssdev
            ss=sstate.(vname);
        end
        iter=iter+1;
        subplot(5,4,iter)
        plot('0:16',myirfs.(shock).(vname)/ss,'linewidth',2)
        title([myvtex{ivar},' to ',myshtex{ishock}])
    end
end
