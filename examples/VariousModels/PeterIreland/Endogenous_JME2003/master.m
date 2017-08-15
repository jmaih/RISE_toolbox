%% housekeeping
clear
close all
clc
%% RISE the model and assign the steady state file
linear=~true;
if linear
    m=rise('emosp_linear','rise_flags',{'original',false},...
        'solve_linear',true);
    % N.B: original =false removes the variables that are both
    % predetermined and forward-looking, making it easier for Peter Ireland
    % to set up the matrices of his system.
else
    m=rise('emosp','steady_state_file','sstate_file',...
        'steady_state_imposed',true);
end
%% Select the model

sticky_price_model=true;
%% get the parameters
[p,priors]=create_parameters(sticky_price_model);
%% push the baseline calibration
m=set(m,'parameters',p);
%% do various things (solving, irfs, vardec, simulation, etc...)
clc
[ms,retcode]=solve(m);
ms.print_solution
%% Select the sample to use for estimation/filtering
choice=2;
samples={
    [],'1979Q2'
    '1979Q3',[]
    [],[]
    };
start_date=samples{choice,1};
end_date=samples{choice,2};
%% load the data
data=create_data(start_date,end_date);

%% estimate the model
log_vars={};
if ~linear
    % The nonlinear model creates problem for the computation of the
    % covariance matrix. Here we take a log-expansion of the state
    % variables
    log_vars=get(m,'endo_list(state)'); % log_vars={'A','E','K','M','V','X','Z'};
end
clc
profile off
profile on
ms=estimate(m,'data',data,...
    'estim_priors',priors,...
    'estim_start_date',start_date,...
    'solve_log_approx_vars',log_vars);%,'kf_tol',0
profile off
profile viewer
%% Redo filtering?
% clc
% [mfilt,LogLik,Incr,retcode]=filter(m,'data',data,'estim_start_date',start_date);%,'kf_tol',0
%% Take a log-linear approximation from a linear approximation
ml=rise('emosp_linear','rise_flags',{'original',true});
mnl=rise('emosp','steady_state_file','sstate_file');
%% Solve and print the solution
clc
M=set([ml,...
    set(mnl,'solve_log_approx_vars',['C',get(mnl,'endo_list(state)')])...
    ],...
    'parameters',p);

M=solve(M);
M.print_solution(['C',M(1).observables.name])%

