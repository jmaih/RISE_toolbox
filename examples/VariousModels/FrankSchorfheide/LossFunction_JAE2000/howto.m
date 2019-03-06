%% housekeeping
clear
clc
close all
%% "rise" the model
m=rise('fs2000_rise');
%% pushing the particular vector (different from the initial estimation point)
% Note we have to declare the regime of the parameter

calibration=struct();
	calibration.alp=0.330;
	calibration.bet=0.990;
	calibration.gam=0.003;
	calibration.mst=1.011;
	calibration.rho=0.700;
	calibration.psi=0.787;
	calibration.del=0.020;
	calibration.sig_a=0.014;
	calibration.sig_m=0.005;
m=set(m,'parameters',calibration); 

%% solve the model and 
m=solve(m);

%% print the solution 

m.print_solution

disp('This should be identical to the dynare solution (see folder dynare_version)')

%% lets do some estimation: we need the data
% get the data from the dynare_version folder

run dynare_version\fsdat_simul

% rise needs data to be passed as a ts object with dates and
% names of the variables. So we proceed to constructing the database

% we use the same start date as in the Schorfheide paper even though the
% data are not the same
startdate='1950q1';
database=ts(startdate,[gp_obs,gy_obs],{'gp_obs','gy_obs'});

%% pass the data to the rise object

% 192 observations are used in dynare
end_date= obs2date(startdate,192);
% the loglinear option of dynare implies:
% 1- we have to take the log of our data. 
% 2- we have to exponentiate the corresponding variables in the rise model
% file. Dynare does this by taking the log of the steady state during
% estimation and this does not ring too transparent to me.

% In rise, we can just take the log of the whole database. This is  
% what we do in passing the data. But we also need to tell ts
% that the logged variables should have the same names as the original
% ones. In order to do that, we add set the flag to true when taking the
% log.
vnames=database.varnames;
database=log(database);
database.varnames=vnames;
m=set(m,'data',database,'estim_end_date',end_date);

%% estimate the model
profile off
profile on
m=estimate(m);%,'optimizer',@csminwellwrap,'debug',true
profile off
profile viewer

%% do posterior simulation
[objective,lb,ub,x0,SIG]=pull_objective(m);

SIG=utils.cov.nearest(SIG);

draws_mcmc = 1000; % number of parameter draws through MCMC.
ndraws_burnin = floor(0.1*draws_mcmc);
mcmc_options=struct('burnin',ndraws_burnin,'N',draws_mcmc,'thin',1,...
    'nchain',2);
Results=mh_sampler(objective,lb,ub,mcmc_options,x0,SIG);

%% update the description of the parameters
m=set(m,'tex_name',...
    {
    'alp','$\alpha$'
    'bet','$\beta$' 
    'gam','$\gamma$' 
    'rho','$\rho$' 
    'psi','$\psi$'
	'del','$\delta$' 
    'sig_a','$\sigma_a$' 
    'sig_m','$\sigma_m$'
    });
%% plot priors, posteriors, priors and posteriors
plot_priors(m)
plot_posteriors(m,Results)
plot_priors_and_posteriors(m,Results)
%% plot priors, posteriors, priors and posteriors for a subset of parameters
close all
myparams={'alp','gam','psi'};
plot_priors(m,myparams)
plot_posteriors(m,Results,myparams)
plot_priors_and_posteriors(m,Results,myparams)
%% check curvature at the mode
mode_curvature(m)

%% check curvature at the mode for a subset of parameters
mode_curvature(m,myparams)

%% check curvature at the mode for a subset of parameters
profile off
profile on
mode_curvature(m,myparams)
profile off
profile viewer
