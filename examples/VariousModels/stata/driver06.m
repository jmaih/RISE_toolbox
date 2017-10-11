%% housekeeping
gentle_clear()
close all
clc

%% load and transform the data
db=xlsread('usmacro2.xlsx');

db=ts('1955Q1',db(:,4:end),{'Y','P','R','C','N','I','E'});

db=pages2struct(db);

%% Rise the madel

m=rise('model06');


%% priors
priors=struct();
priors.beta={0.98,0.1,0.999};
priors.kappa={0.1,0.001,5};
priors.rhoz={0.5,-0.999,0.999};
priors.rhor={0.5,-0.999,0.999};
priors.rhou={0.5,-0.999,0.999};
priors.sigu={0.1,0.001,10};
priors.sigz={0.1,0.001,10};

%% estimate model

mest=estimate(m,'data_demean',true,'data',db,'priors',priors,...
    'estim_start_date','1955Q1','estim_end_date','2015Q4');

