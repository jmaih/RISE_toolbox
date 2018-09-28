%% housekeeping
gentle_clear()
close all
clc

%% load and transform the data
db=xlsread('usmacro2.xlsx');

db=ts('1955Q1',db(:,4:end),{'Y','P','R','C','N','I','E'});

db=pages2struct(db);

%% Rise the madel

m=rise('model02');

%% fixed parameters

m=set(m,'parameters',{'beta',0.96});

%% Estimate the madel

priors=struct();

priors.kappa={0.2,0.001,2};
priors.psi={1.5,1,3};
priors.rhou={0.5,0,0.999};
priors.rhog={0.5,0,0.999};
priors.sigu={2,0.0001,5}; 
priors.sigg={0.5,0.0001,5};

[mest,filtration]=estimate(m,'priors',priors,'data_demean',true,'data',db,...
    'estim_start_date','1955Q1','estim_end_date','2015Q4');

%% print solution at the mode

print_solution(mest)

%% Plot Inflation vs one-step ahead

fs=filtration.smoothed_variables;

ff=filtration.filtered_variables;

figure('name','Inflation vs predicted (one-step) inflation')
% The mean is added back to the observables after filtering
plot([fs.P,ff.P])


