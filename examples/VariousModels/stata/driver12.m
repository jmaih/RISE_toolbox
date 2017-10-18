% A model with two unidentified parameters

%% housekeeping
gentle_clear()
close all
clc

%% load and transform the data
db=xlsread('usmacro2.xlsx');

db=ts('1955Q1',db(:,4:end),{'Y','P','R','C','N','I','E'});

db=pages2struct(db);

%% Rise the madel

m=rise('model12');

%% Set gam to 1

m=set(m,'parameters',{'gam',1});

%% priors
priors=struct();
priors.beta={0.98,0.1,0.999};
priors.kappa={0.1,0.001,5};
priors.gam={0.5,0.001,5};
priors.rhoz={0.5,-0.999,0.999};
priors.rhou={0.5,-0.999,0.999};
priors.sigu={0.1,0.001,15};
priors.sigz={0.1,0.001,15};

%% estimate model fixing gam
mest1=estimate(m,'data_demean',true,'data',db,...
    'priors',rmfield(priors,'gam'),...
    'estim_start_date','1955Q1','estim_end_date','2015Q4');

%% investigate identification

J1=dsge_tools.identification(mest1,[],[],5);

size(J1,2)-rank(full(J1),1e-6)

%% estimate model without fixing gam
mest2=estimate(m,'data_demean',true,'data',db,...
    'priors',priors,...
    'estim_start_date','1955Q1','estim_end_date','2015Q4');

%% investigate identification

J2=dsge_tools.identification(mest2,[],[],5);

size(J2,2)-rank(full(J2),1e-6)


