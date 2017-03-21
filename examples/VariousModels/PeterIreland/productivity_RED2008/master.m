%% housekeeping
clc

%% rise the model

m=rise('pdk','rise_flags',{'Sektors',{'C','I'}},...
    'steady_state_file','sstate_model',...
    'saveas',true);

%% get the parameters

[p,priors]=create_parameters(~true);

m=set(m,'parameters',p);

%% get the data

data=create_data();

%% estimate model
clc

ms=estimate(m,'data',data,'estim_priors',priors,'optimizer','fminunc');
