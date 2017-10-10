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
%% Choose a model type: see cell "create the structural VAR model" below
model_type=4;

%% set up the restrictions
close()

% create restrictions on parameters as well as markov chains
%------------------------------------------------------------
switch model_type
    case 0 
        % constant-parameter model
        [lin_restr,nonlin_restr,markov_chains]=create_restrictions_and_markov_chains0();
    case 1 
        % Coefficients are switching regimes across all equations
        % (synchronized case) 
        [lin_restr,nonlin_restr,markov_chains]=create_restrictions_and_markov_chains1();
    case 2 
        % Coefficients and variances have different chains, different
        % regimes, and different durations 
        [lin_restr,nonlin_restr,markov_chains]=create_restrictions_and_markov_chains2();
    case 3 
        % Only coefficients in monetary policy equation are changing
        [lin_restr,nonlin_restr,markov_chains]=create_restrictions_and_markov_chains3();
    case 4 
        % Only variance in monetary policy equation is changing
        [lin_restr,nonlin_restr,markov_chains]=create_restrictions_and_markov_chains4();
    case 5 
        % Both coefficients and variances in monetary policy equation
        % change with two independent Markov processes 
        [lin_restr,nonlin_restr,markov_chains]=create_restrictions_and_markov_chains5();
    case 6 % ok
        % Only variances in ALL three equations switch
        [lin_restr,nonlin_restr,markov_chains]=create_restrictions_and_markov_chains6();
    otherwise
        error('the coded model types are 0, 1, 2, 3, 4 and 6')
end

%% Create the VAR

nlags=2;

exog={};

panel=[];

constant=true;

% first we create a template structure
% ------------------------------------
sv0=svar(varlist,exog,nlags,constant,panel,markov_chains);

%% set priors

prior=svar.prior_template();

prior.type='sz';

% prior.L1=0.1/2;
% 
% prior.coefprior=0.5;

is_prior=~true;

%% Find posterior mode
clc

sv=sv0;

if is_prior
    
    sv=set(sv,'prior',prior);
    
end

sv=estimate(sv,{'1960Q1','2015Q2'},'data',db,...
    'linear_restrictions',[lin_restr;nonlin_restr]);

%% Printing estimates

print_structural_form(sv)

%% Printing solution
clc

print_solution(sv)

%% plot smoothed state and regime probabilities

close all

plot_probabilities(sv)

%% plots probabilities against data
close all

plot_data_against_probabilities(sv,'state')

%% Impulse responses

myirfs=irf(sv);

%% Posterior sampling

%% Marginal data density

%% Out-of sample forecasts

%% Conditional forecast


