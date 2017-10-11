%% housekeeping
gentle_clear()
close all
clc

%% load and transform the data
db=xlsread('usmacro2.xlsx');

db=ts('1955Q1',db(:,4:end),{'Y','P','R','C','N','I','E'});

db=pages2struct(db);

%% Rise the madel

m=rise('model03');

%% fixed parameters

p={'beta',0.96
    'eta',1
    'alpha',0.3
    'delta',0.025
    'phi1',0.2
    'phi2',0.6
    'rhoz',0.8
    'rhog',0.3
    'sigz',1
    'sigg',1};


m=set(m,'parameters',p);

%% solve

m=solve(m);

%% compute likelihood
[mfilt,LogLik,Incr,retcode]=filter(m,'data_demean',true,'data',db,...
    'estim_start_date','1955Q1','estim_end_date','2015Q4');

fprintf(1,'%0.4f\n',LogLik)

%% print solution at the mode

print_solution(m)

%% Impulse responses

myirfs=irf(m);

quick_irfs(m,myirfs,[],[],[3,3],'0:8')


