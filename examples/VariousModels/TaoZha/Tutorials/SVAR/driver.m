%% housekeeping
close all
home()
%% Choose a model type: see cell "create the structural VAR model" below
model_type=2;

%% create dataset
do_plot=true;

[db,varlist]=create_dataset(do_plot);

%% create the structural VAR model
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
        [restrictions,tpl]=create_restrictions_and_markov_chains0(tpl);
    case 1
        % Coefficients are switching regimes across all equations
        % (synchronized case) 
        [restrictions,tpl]=create_restrictions_and_markov_chains1(tpl);
    case 2
        % Coefficients and variances have different chains, different
        % regimes, and different durations 
        [restrictions,tpl]=create_restrictions_and_markov_chains2(tpl);
    case 3
        % Only coefficients in monetary policy equation are changing
        [restrictions,tpl]=create_restrictions_and_markov_chains3(tpl);
    case 4
        % Only variance in monetary policy equation is changing
        [restrictions,tpl]=create_restrictions_and_markov_chains4(tpl);
    case 5
        % Both coefficients and variances in monetary policy equation
        % change with two independent Markov processes 
        [restrictions,tpl]=create_restrictions_and_markov_chains5(tpl);
    otherwise
        error('the coded model types are 0, 1, 2, 3, 4 and 5')
end

% finally we create a svar object by pushing the structure into svar
%--------------------------------------------------------------------
m=svar(tpl,'data',db,'estim_linear_restrictions',restrictions);

%% find posterior mode
profile off
profile on
mest=estimate(m,'estim_start_date','1960Q1');
profile off
profile viewer

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
clc
mdd_algorithms={'bridge','mhm','mueller','swz','is','ris','cj'};
nalgos=numel(mdd_algorithms);
log_mdd=zeros(1,nalgos);
extras=cell(1,nalgos);
howlong=zeros(1,nalgos);
for ichoice=1:nalgos
    tic
    [log_mdd(ichoice),extras{ichoice}] = mcmc_mdd(Results.pop,lb,ub,...
        struct('log_post_kern',objective,...
        'algorithm',mdd_algorithms{ichoice},...
        'L',500,'debug',true));
    howlong(ichoice)=toc;
    fprintf(1,'log MDD(%s): %0.4f time: %0.4f seconds\n',...
        mdd_algorithms{ichoice},log_mdd(ichoice),...
        howlong(ichoice));
end

%% Impulse responses
myirfs=irf(mest);

for ishock=1:numel(mest.exogenous.name)
    shock_name=mest.exogenous.name{ishock};
    figure('name',['IRF to a ',shock_name,' shock'])
    for ivar=1:numel(varlist)
        vname=varlist{ivar};
        subplot(3,1,ivar)
        plot(myirfs.(shock_name).(vname),'linewidth',2)
        title(vname)
    end
end
%% Out-of sample forecasts

mycast=forecast(mest);

%% Conditional forecast on ygap: parameter uncertainty only...
ndraws=200;
[fkst,bands]=do_conditional_forecasts(mest,db,Results.pop,ndraws);

%%
figure('name','Conditional forecasts on ygap');
for ivar=1:mest.endogenous.number
    subplot(mest.endogenous.number,1,ivar)
    vname=mest.endogenous.name{ivar};
    plot(bands.(vname),'linewidth',2)
    title(vname)
    if ivar==3
        leg=legend(get(bands.(vname),'varnames'));
        set(leg,'interp','none','Orientation','horizontal',...
            'location','SouthOutside')
    end
end
