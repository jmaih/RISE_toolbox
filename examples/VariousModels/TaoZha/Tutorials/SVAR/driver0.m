%% --------------------- constant-parameter model --------------------- %%

%% housekeeping
clear 
close all
clc()

%% Create dataset
clc

do_plot=true;

scale=100;

[db,varlist0]=create_dataset(scale,do_plot);

varlist=fieldnames(varlist0);

%% Create the VAR

clc

nlags=2;

exog={};

constant=true;

sv0=svar(varlist,exog,nlags,constant);

%% set up restrictions

% syntax is alag(eqtn,vname)
%-------------------------------
lin_restr={
    % first equation or "FFR" equation
    %----------------------------------
    'a1(1,pi)=0'
    'a2(1,pi)=0'
    'a1(1,ygap)=0'
    'a2(1,ygap)=0'
    'a2(1,FFR)=0'
    % second equation or "pi" equation
    %----------------------------------
    'a0(2,FFR)=0'
    'a1(2,FFR)=0'
    'a2(2,FFR)=0'
    'a1(2,ygap)=0'
    'a2(2,ygap)=0'
    % third equation or "ygap" equation
    %-----------------------------------
    'a1(3,FFR)=0'
    'a2(3,FFR)=0'
    'a1(3,pi)=0'
    'a2(3,pi)=0'
    'a0(3,pi)+a0(3,FFR)=0'
    };
nonlin_restr={
    'a0(3,FFR)>=0'
    'a1(1,FFR)>=0'
    'a1(1,FFR)<=1'
    };

restrictions=[lin_restr;nonlin_restr];

%% set priors 

var_prior=svar.prior_template();

var_prior.type='sz';

prior=struct('var',var_prior);


%% Find posterior mode
clc

sv=sv0;

sv=estimate(sv,db,{'1960Q1','2015Q2'},prior,restrictions);

%% estimates
clc

pmode=posterior_mode(sv)

%% Printing estimates
clc

print_structural_form(sv)

%% Printing solution
clc

print_solution(sv)

%% plot smoothed state and regime probabilities
clc
close all

plot_probabilities(sv)

%% plots probabilities against data
close all

plot_data_against_probabilities(sv,'regime')

%% Impulse responses

myirfs=irf(sv);

%% Posterior sampling

%% Marginal data density

%% Out-of sample forecasts

%% Conditional forecast


