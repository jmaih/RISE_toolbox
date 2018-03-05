%% --------------------- all coefficients switch --------------------- %%
%
% Coefficients are switching regimes across all equations
% (synchronized case)

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

%% set up Markov chains
% We add a Markov chain to the template
%----------------------------------------
% N.B: The chain controls the coefficients of all equations but not the
% variance, which remains constant.
markov_chains=struct('name','syncoef',...
    'number_of_states',2,...
    'controlled_parameters',{{'c','a0','a1','a2'}},...
    'endogenous_probabilities',[],...
    'probability_parameters',[]);

%% Create the VAR

clc

nlags=2;

exog={};

constant=true;

panel=[];

sv0=svar(varlist,exog,nlags,constant,panel,markov_chains);

%% set up restrictions

% syntax is alag(eqtn,vname,chain_name,state)
%------------------------------------------------
lin_restr=cell(0,1);
nonlin_restr=cell(0,1);

numberOfStates=markov_chains.number_of_states;

for istate=1:numberOfStates
    
    mystate=int2str(istate);
    
    lin_restr=[lin_restr
        {
        % first equation or "FFR" equation:
        %----------------------------------
        ['a1(1,pi,syncoef,',mystate,')=0']
        ['a2(1,pi,syncoef,',mystate,')=0']
        ['a1(1,ygap,syncoef,',mystate,')=0']
        ['a2(1,ygap,syncoef,',mystate,')=0']
        ['a2(1,FFR,syncoef,',mystate,')=0']
        % second equation or "pi" equation
        %----------------------------------
        ['a0(2,FFR,syncoef,',mystate,')=0']
        ['a1(2,FFR,syncoef,',mystate,')=0']
        ['a2(2,FFR,syncoef,',mystate,')=0']
        ['a1(2,ygap,syncoef,',mystate,')=0']
        ['a2(2,ygap,syncoef,',mystate,')=0']
        % third equation or "ygap" equation
        %-----------------------------------
        ['a1(3,FFR,syncoef,',mystate,')=0']
        ['a2(3,FFR,syncoef,',mystate,')=0']
        ['a1(3,pi,syncoef,',mystate,')=0']
        ['a2(3,pi,syncoef,',mystate,')=0']
        ['a0(3,pi,syncoef,',mystate,')+a0(3,FFR,syncoef,',mystate,')=0']
        }
        ]; %#ok<AGROW>
    nonlin_restr=[nonlin_restr
        {
        ['a0(3,FFR,syncoef,',mystate,')>=0']
        ['a1(1,FFR,syncoef,',mystate,')>=0']
        ['a1(1,FFR,syncoef,',mystate,')<=1']
        }
        ]; %#ok<AGROW>
end

restrictions=[lin_restr;nonlin_restr];

%% set priors 

var_prior=svar.prior_template();

var_prior.type='sz';

switch_prior=struct();
switch_prior.syncoef_tp_1_2={0.5,0.1,0.3,'beta'};
switch_prior.syncoef_tp_2_1={0.5,0.1,0.3,'beta'};

prior=struct();

prior.var=var_prior;

prior.nonvar=switch_prior;

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


