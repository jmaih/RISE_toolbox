%% housekeeping
clc
%% RISE the model
linear=~true;
if linear
    m=rise('targets_lin');
else
    m=rise('targets','steady_state_file','sstate_model');
end
%% get the parameters
version='a';
% 'a' % model with an endogenous inflation target.
% 'b' % model with an exogenous inflation target.
% 'd' % model with backward-looking price setting.
start_at_mode=false;

[p,priors]=create_parameters(version,linear,start_at_mode);

%% push the parameters
m=set(m,'parameters',p);

%% create data
[data]=create_data();

%% estimate the model
ms=estimate(m,'data',data,'estim_priors',priors);

%% Impulse responses
myirfs=irf(ms);

%% plot responses
myvars={'Y','PAI','R','X'};
locs=locate_variables(myvars,ms.endogenous.name);
myvtex=ms.endogenous.tex_name(locs);
sstate=get(ms,'sstate');
ssdev=false;

shocks={'EPS_A','EPS_E','EPS_Z','EPS_V','EPS_PAI'};%m.exogenous.name;
locs=locate_variables(shocks,ms.exogenous.name);
myshtex=ms.exogenous.tex_name(locs);
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
