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

m=rise('trinity');

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
parameter_names={m.estimation.priors.name};
lb=vertcat(m.estimation.priors.lower_bound);
ub=vertcat(m.estimation.priors.upper_bound);

%%
nsim_or_draws=1000; %price_puzzle

check_behavior=@(x)is_has_solution(m,parameter_names,x);

procedure='latin_hypercube';

obj=mcf(check_behavior,nsim_or_draws,lb,ub,parameter_names,procedure);

%% smirnov test for equality of distributions

% hfig=utils.plot.multiple(@(x)cdf_plot(obj,x),parameter_names,...
%     'Smirnov test of equality of distributions',3,3)
% 
cdf_plot(obj)

%% correlation patterns in the behavior sample

correlation_patterns_plot(obj,[],'behave')

%% correlation patterns in the behavior sample

correlation_patterns_plot(obj,[],'non-behave')

%% scatter plot of the significant correlations

scatter(obj,[],'non-behave')
