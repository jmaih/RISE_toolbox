%% housekeeping
gentle_clear()
close all
clc

%% load and transform the data
db=xlsread('usmacro2.xlsx');

db=ts('1955Q1',db(:,4:end),{'Y','P','R','C','N','I','E'});

db=pages2struct(db);

%% Rise the madel

m=rise('model11');

%% priors
priors=struct();
priors.alpha={0.1,-1,1};
priors.rhog={0.5,0,0.999};
priors.rhogz={0.5,-0.999,0.999};
priors.rhoz={0.5,0,0.999};
priors.sigg={0.1,0.001,15};
priors.sigz={0.1,0.001,15};
%% estimate model

% It is a mystery to me why MFI works here and not Klein...

mest=estimate(m,'data_demean',true,'data',db,'priors',priors,...
    'estim_start_date','1955Q1','estim_end_date','2015Q4',...
    'solver','mfi');

