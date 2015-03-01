%% housekeeping
clear all
close all
clc
%% RISE the model and assign the steady state file
linear=true;
if linear
    m=rise('ireland2003_linear');
else
    m=rise('ireland2003','steady_state_file','ireland2003_sstate');
end
%% get the parameters
[p,priors]=ireland2003_parameterization();
%% push the baseline calibration
m=set(m,'parameters',p);
%% do various things (solving, irfs, vardec, simulation, etc...)
clc
[ms,retcode]=solve(m);
ms.print_solution
%% load the data
data=ireland2003_create_data();
%%
profile off
profile on
ms=estimate(m,'data',data,'kf_tol',0,'estim_priors',priors,'estim_start_date','1979Q3');
profile off
profile viewer
%% Redo filtering?
clc
[mfilt,LogLik,Incr,retcode]=filter(m,'data',data,'kf_tol',0,'estim_end_date','1979Q2');