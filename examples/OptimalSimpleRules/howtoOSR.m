%% housekeeping
close all
clear
clc

%% read the model 

m=rise('Canonical_osr');

%% Baseline parameters
p={
    'beta_lag' 	 ,0.5   
	'beta_lead'	 ,0.4	 
	'beta_r'  	 ,0.9
	'lamb_lag'	 ,0.8 
	'lamb_lead'    ,0.1 
	'lamb_y'  	 ,0.3 
	'sigi'   		 ,0.5	 
	'sigpai' 		 ,0.5	 
	'sigy'   		 ,0.5
	'gam_lag' 	 ,0.6	 
	'gam_y'   	 ,0.5
    };
m=set(m,'parameters',p); 

%% priors
priors=struct();
priors.gam_lag={0.6000, 0.0000, 1};	 
priors.gam_y={0.5000, 0.0000, 10};	

%% Loss function
% no need to specify a list of observable variables

% specify the planner's objective. 
Loss='-.5*(1*PAI^2+.3*Y^2+0.9*DI^2)';%{discount = 0.99}
% The discount factor does not play a role. So I can... neglect it.
%planner_objective{discount = 0.99} -.5*(1*PAI^2+.3*Y^2+0.9*DI^2);	

%% Posterior maximization
simuls=[];

[mest]=optimal_simple_rule(m,Loss,simuls,'priors',priors)

%% estimate with fmincon
m=osr(m);

