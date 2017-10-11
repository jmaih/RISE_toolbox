%% housekeeping
gentle_clear()
close all
clc

%% load and transform the data
db=xlsread('usmacro2.xlsx');

db=ts('1955Q1',db(:,4:end),{'Y','P','R','C','N','I','E'});

db=pages2struct(db);

%% Rise the madel

m=rise('model10');

%% priors
priors=struct();
priors.beta={0.98,0.1,0.999};
priors.kappa={0.1,0.001,5};
priors.psi={0.5,0.001,1};
priors.rhou={0.5,-0.999,0.999};
priors.rhog={0.5,-0.999,0.999};
priors.rhoe={0.5,-0.999,0.999};
priors.sigu={0.1,0.001,15};
priors.sigg={0.1,0.001,15};
priors.sige={0.1,0.001,15};
%% estimate model

mest=estimate(m,'data_demean',true,'data',db,'priors',priors,...
    'estim_start_date','1973Q2','estim_end_date','2015Q4');

