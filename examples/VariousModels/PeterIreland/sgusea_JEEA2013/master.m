%% Housekeeping
clc

%% RISE the model

m=rise('sgusea1',...
    'rise_flags',struct('CountryNames',{{'H','F'}}),...
    'rise_save_macro',true);

%% Assign steady state file

m=set(m,'steady_state_file','sstate_model');

%% create parameters

[p,priors]=create_parameters(true);

m=set(m,'parameters',p);

%% Bring in the data

data=create_data();

%% Estibrate the model
clc
ms=estimate(m,'data',data,'estim_priors',priors);
%%
[mfilt,loglik]=filter(m,'data',data);

