%% housekeeping
clear
close all
clc

%% load and transform the data
M=importdata('FEDFUNDS.csv');

Q=importdata('GDPDEF.csv');

db=struct();

db.GDPDEF=ts('1947Q1',Q.data);

db.FEDFUNDS=aggregate(ts('1954M7',M.data),'Q');

db.P=400*log(db.GDPDEF/db.GDPDEF{-1});

db.R=db.FEDFUNDS;

%% 

m=rise('model1');

%% calibration and irfs
p=struct();
p.beta=0.5112881;
p.kappa=0.1696296;
p.rhou=0.6989189;
p.rhog=.9556407;
p.sigu=2.317589;
p.sigg=.6147348;

myirfs=irf(m,'parameters',p);

%%

quick_irfs(m,myirfs,[],[],[3,2])

%%

priors=struct();

priors.beta={0.96,0.1,0.999};
priors.kappa={0.2,0.001,2};
priors.rhou={0.5,0,0.999};
priors.rhog={0.5,0,0.999};
priors.sigu={2,0.0001,5}; 
priors.sigg={0.5,0.0001,5};

mest=estimate(m,'priors',priors,'data_demean',true,'data',db,...
    'estim_start_date','1954Q3','estim_end_date','2016Q4');
