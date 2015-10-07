%% housekeeping
close all
home()
%% Choose a model type: see cell "create the structural VAR model" below
model_type=5;

%% Create dataset
do_plot=true;
scale=100;

[db,varlist]=create_dataset(scale,do_plot);

%% Create the structural VAR model
close()

% first we create a template structure
% ------------------------------------
tpl=svar.template();

% we update the fields of the structure
% --------------------------------------
tpl.endogenous=varlist;
tpl.nlags=2;

% create restrictions on parameters as well as markov chains
%------------------------------------------------------------
switch model_type
    case 0
        % constant-parameter model
        [lin_restr,nonlin_restr,tpl]=create_restrictions_and_markov_chains0(tpl);
    case 1
        % Coefficients are switching regimes across all equations
        % (synchronized case) 
        [lin_restr,nonlin_restr,tpl]=create_restrictions_and_markov_chains1(tpl);
    case 2
        % Coefficients and variances have different chains, different
        % regimes, and different durations 
        [lin_restr,nonlin_restr,tpl]=create_restrictions_and_markov_chains2(tpl);
    case 3
        % Only coefficients in monetary policy equation are changing
        [lin_restr,nonlin_restr,tpl]=create_restrictions_and_markov_chains3(tpl);
    case 4
        % Only variance in monetary policy equation is changing
        [lin_restr,nonlin_restr,tpl]=create_restrictions_and_markov_chains4(tpl);
    case 5
        % Both coefficients and variances in monetary policy equation
        % change with two independent Markov processes 
        [lin_restr,nonlin_restr,tpl]=create_restrictions_and_markov_chains5(tpl);
    case 6
        % Only variances in ALL three equations switch
        [lin_restr,nonlin_restr,tpl]=create_restrictions_and_markov_chains6(tpl);
    otherwise
        error('the coded model types are 0, 1, 2, 3, 4 and 6')
end

% finally we create a svar object by pushing the structure into svar
%--------------------------------------------------------------------
m=svar(tpl,'data',db,'estim_linear_restrictions',lin_restr,...
    'estim_nonlinear_restrictions',nonlin_restr);

%% Find posterior mode

mest=estimate(m,'estim_start_date','1960Q1');

%% Printing out A0hat, A1hat, and A2hat in a form compatible with Zha''s original Matlab code

output_in_original_zha_matlab_code(mest,scale)

% uncomment this if you want to capture output
% [A0hat,A1hat,A2hat]=output_in_original_zha_matlab_code(mest,scale)

%% Markov chain Monte Carlo
% Note that because of the linear restrictions, not all parameters are
% estimated. Hence, the effective number of estimated parameters is smaller
% than the number of parameters declared by the user. This is reflected in
% the dimensions of lb,ub,x0 and SIG below. The user does not have to be
% concerned about those.
[objective,lb,ub,x0,SIG]=pull_objective(mest);

ndraws_mcmc         = 1500;  % number of parameter draws through MCMC.
ndraws_burnin       = floor(0.1*ndraws_mcmc); % number of parameter draws to be burned
mcmc_options=struct('burnin',ndraws_burnin,'N',ndraws_mcmc,'thin',1);

Results=mh_sampler(objective,lb,ub,mcmc_options,x0,SIG);

%% Marginal data density
% pick yours: 'bridge','mhm','mueller','swz','is','ris','cj'
tic
log_mdd = mcmc_mdd(Results.pop,lb,ub,...
    struct('log_post_kern',objective,... % function to MINIMIZE !!!
    'algorithm','swz',... % MDD algorithm
    'L',2000 ... % Number of i.i.d. draws of the proposal density function
));
minutes_MDD_Took = toc/60;

laplace_approx=mest.estimation.posterior_maximization.log_marginal_data_density_laplace;
disp(['log(MDD): Laplace Approximation: ',num2str(laplace_approx)])
disp(['log(MDD): Sims-Waggoner-Zha Approximation: ',num2str(log_mdd)])
disp([' Minutes Sims-Waggoner-Zha Approximation took: ',num2str(minutes_MDD_Took)])

%% Impulse responses
myirfs=irf(mest);

do_plot_irfs(myirfs,mest);

%% do smoothed probabilities

do_plot_smoothed_probabilities(mest)

%% Out-of sample forecasts at the mode

mycast=forecast(mest);

do_plot_unconditional_forecasts(mycast,mest)

%% Conditional forecast on ygap
% conditional information
%-------------------------
ygap=scale*(-0.025:0.005:8*0.005).';
cond_db=struct('ygap',ts('2015Q3',ygap));

% options for the exercise
%--------------------------
myoptions=struct('cbands',[10,20,50,80,90],'do_plot',true,'nsteps',20,...
    'param_uncertainty',false,'shock_uncertainty',true,'ndraws',200,...
    'simul_regime',[]);

% do it
%-------
tic
[fkst,bands,hdl]=do_conditional_forecasts(mest,db,cond_db,Results.pop,myoptions);
fprintf('\n\n Computing conditional forecasts took %0.4f minutes\n\n',toc/60);
%% Median conditional forecasts
% we take advantage of the fact that in the bands above we specified 50 in
% the bands above 
%--------------------------------------------------------------------------
figure('name','Median Conditional Forecasts');
nvars=mest.endogenous.number;
for ivar=1:nvars
    thisname=mest.endogenous.name{ivar};
    subplot(nvars,1,ivar)
    plot(bands.(thisname)('ci_50'),'linewidth',2)
    title(mest.endogenous.tex_name{ivar})
end


