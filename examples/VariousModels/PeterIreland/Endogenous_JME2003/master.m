%% housekeeping
clear all
close all
clc
%% RISE the model and assign the steady state file
linear=true;
if linear
    m=rise('emosp_linear','rise_flags',{'original',false});
else
    m=rise('emosp','steady_state_file','sstate_file');
end
%% get the parameters
[p,priors]=create_parameters();
%% push the baseline calibration
m=set(m,'parameters',p);
%% do various things (solving, irfs, vardec, simulation, etc...)
clc
[ms,retcode]=solve(m);
ms.print_solution
%% load the data
data=create_data();
%%
profile off
profile on
ms=estimate(m,'data',data,'kf_tol',0,'estim_priors',priors,'estim_start_date','1979Q3');
profile off
profile viewer
%% Redo filtering?
clc
[mfilt,LogLik,Incr,retcode]=filter(m,'data',data,'kf_tol',0,'estim_end_date','1979Q2');