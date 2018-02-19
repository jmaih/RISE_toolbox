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

model_type=0;

%% set up the restrictions
close()

% create restrictions on parameters as well as markov chains
%------------------------------------------------------------
switch model_type
    case 0 
        % constant-parameter model
        [restrictions,markov_chains,switch_prior]=create_restrictions_and_markov_chains0();
    case 1 
        % Coefficients are switching regimes across all equations
        % (synchronized case) 
        [restrictions,markov_chains,switch_prior]=create_restrictions_and_markov_chains1();
    case 2 
        % Coefficients and variances have different chains, different
        % regimes, and different durations 
        [restrictions,markov_chains,switch_prior]=create_restrictions_and_markov_chains2();
    case 3 
        % Only coefficients in monetary policy equation are changing
        [restrictions,markov_chains,switch_prior]=create_restrictions_and_markov_chains3();
    case 4 
        % Only variance in monetary policy equation is changing
        [restrictions,markov_chains,switch_prior]=create_restrictions_and_markov_chains4();
    case 5 
        % Both coefficients and variances in monetary policy equation
        % change with two independent Markov processes 
        [restrictions,markov_chains,switch_prior]=create_restrictions_and_markov_chains5();
    case 6 % ok
        % Only variances in ALL three equations switch
        [restrictions,markov_chains,switch_prior]=create_restrictions_and_markov_chains6();
    otherwise
        error('the coded model types are 0, 1, 2, 3, 4 and 6')
end

%% Create the VAR

clc

nlags=2;

exog={};

panel=[];

constant=true;

% first we create a template structure
% ------------------------------------
sv0=svar(varlist,exog,nlags,constant,panel,markov_chains);

%% set priors % prior=[];

var_prior=svar.prior_template();

var_prior.type='sz';

prior=struct('var',var_prior,'nonvar',switch_prior);

is_prior=true;

if ~is_prior
    
    prior=rmfield(prior,'var');
    
end

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


