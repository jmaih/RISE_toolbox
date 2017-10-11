%% housekeeping
gentle_clear()
close all
clc

%% load and transform the data
db=xlsread('rates2.xlsx');

db=ts('1947Q1',db(:,1:2),{'gdpdef','R'});

db=pages2struct(db);

db.P=400*log(db.gdpdef/db.gdpdef{-1});

%% Rise the madel

m=rise('model01');

%% Estimate the madel

priors=struct();

priors.beta={0.96,0.1,0.999};
priors.kappa={0.2,0.001,2};
priors.rhou={0.5,0,0.999};
priors.rhog={0.5,0,0.999};
priors.sigu={2,0.0001,5}; 
priors.sigg={0.5,0.0001,5};

mest=estimate(m,'priors',priors,'data_demean',true,'data',db,...
    'estim_start_date','1954Q3','estim_end_date','2016Q4');

%% print solution at the mode

print_solution(mest)

%% irfs at the mode

% myirfs=irf(m,'parameters',get(mest,'parameters'));

myirfs=irf(mest);

%%

quick_irfs(m,myirfs,[],[],[3,2])
