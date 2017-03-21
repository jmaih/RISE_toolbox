%% Housekeeping
% In this version, one can choose to run model 1 or model 2 by selecting
% the appropriate flag "ncp_shocks". If ncp_shocks==false, we have model 1
% otherwise we have model 2.
% Model 1: nonstationary and cointegrated
% - Neutral technology shocks
% - Investment-specific technology shocks
% Model 2: nonstationary and cointegrated
% - Neutral technology shocks
% - Investment-specific technology shocks
% - Preference shocks
% cointegrated
clc

%% RISE the model
ncp_shocks=true;
m=rise('sgusea12',...
    'rise_flags',struct('CountryNames',{{'H','F'}},'ncp_shocks',ncp_shocks),...
    'saveas',true);

%% Assign steady state file

m=set(m,'steady_state_file','sstate_model12b');

%% create parameters

[p,priors]=create_parameters12(ncp_shocks);

m=set(m,'parameters',p);

%% Bring in the data

data=create_data();

%% Estibrate the model
clc
ms=estimate(m,'data',data,'estim_priors',priors);
%%
[mfilt,loglik]=filter(m,'data',data);

