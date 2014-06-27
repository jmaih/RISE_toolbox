%% housekeeping
close all
clear
clc
%% a few settings

% Only correlations with pval above pattern_cutoff will be considered for
% ploting of the correlation patterns
pattern_cutoff=.1;

% Only correlations with pval above scatter_cutoff will be considered for
% ploting in the scatter plots
scatter_cutoff=.01;

%% rise the model

tnt=rise('trinity');

%% we need a function describing the behavior
% the function will take as first argument, the rise object and possibly
% other arguments
% the main argument returned by the function shall be a logical scalar,
% assuming the value true if the parameter vector is consistent with the
% behavior or false otherwise.
% QUESTION: What is the region of the parameter space which is such that 
% prices increase after a monetary policy shock?
% ANSWER: pass the function price_puzzle to tnt and find the answer.
% But before doing that, let's print the function on screen
clc
type price_puzzle
%%
% the parameters that are consistent with the price puzzle behavior
tnt=set_options(tnt,'monte_carlo_behavioral_function',@price_puzzle);
%% Now we can go ahead and run the Monte Carlo Filtering
[xparam,retcode,could_solve,behave]=monte_carlo_filtering(tnt);
%% correlation patterns in the behavior sample
parameter_names={tnt.estimation.priors.name};
lb=vertcat(tnt.estimation.priors.lower_bound);
ub=vertcat(tnt.estimation.priors.upper_bound);

hh=mcf.plot_correlation_patterns(xparam(:,could_solve),...
    parameter_names,pattern_cutoff,'behavior');

%% correlation patterns in the non-behavior sample
hh2=mcf.plot_correlation_patterns(xparam(:,~could_solve),...
    parameter_names,pattern_cutoff,'non-behavior');
%% Smirnov test for equality of distributions 
mcf.plot_smirnov(xparam,could_solve,lb,ub,parameter_names,4,3)
%% scatter plot of the significant correlations
mcf.plot_correlation_scatter(xparam,could_solve,parameter_names,...
    scatter_cutoff,4,3)
