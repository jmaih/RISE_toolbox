%% housekeeping
gentle_clear()
close all
clc

%% load and transform the data
db=xlsread('usmacro2.xlsx');

db=ts('1955Q1',db(:,4:end),{'Y','P','R','C','N','I','E'});

db=pages2struct(db);

%% Rise the madel

m=rise('model04');

%% fixed parameters

p={'beta',0.96
    'psi',2
    'chi',0.8};

m=set(m,'parameters',p);

%% priors
priors=struct();
priors.kappa={0.1,0.0001,2};
priors.chi={0.1,0.001,1};
priors.psi={1.5,1,10};
priors.rhoe={0.5,0,0.999};
priors.rhou={0.5,0,0.999};
priors.rhog={0.5,0,0.999};
priors.sige={0.1,0.001,10};
priors.sigu={0.1,0.001,10};
priors.sigg={0.1,0.001,10};

%% compute likelihood
mest=estimate(m,'data_demean',true,'data',db,'priors',priors,...
    'estim_start_date','1955Q1','estim_end_date','2015Q4');


%% print solution at the mode

print_solution(mest)

%% Impulse responses

myirfs=irf(mest);

quick_irfs(m,myirfs,[],[],[3,3],'0:8')


