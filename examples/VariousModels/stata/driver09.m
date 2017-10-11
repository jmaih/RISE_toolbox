%% housekeeping
gentle_clear()
close all
clc

%% load and transform the data
db=xlsread('usmacro2.xlsx');

db=ts('1955Q1',db(:,4:end),{'Y','P','R','C','N','I','E'});

db=pages2struct(db);

%% Rise the madel

m=rise('model09');

%% priors
priors=struct();
priors.b1={0.5,0.1,1};
priors.gam={0.1,0.001,5};
priors.h={0.5,0.1,1};
priors.rhow={0.5,-0.999,0.999};
priors.sigw={0.1,0.001,10};
priors.sigr={0.1,0.001,10};

%% estimate model

mest=estimate(m,'data_demean',true,'data',db,'priors',priors,...
    'estim_start_date','1955Q1','estim_end_date','2015Q4');

