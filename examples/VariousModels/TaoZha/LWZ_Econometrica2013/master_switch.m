%% housekeeping
close all
clc
clear

%% add a path

addpath('routines')

%% rise the model

m=rise('lwz_switch_model','solve_linear',true);

%% get the parameters

[p,priors]=create_parameters_switch(true);

m=set(m,'parameters',p);

%% bring in the data

data=create_data();

data=pages2struct(data);

vnames=fieldnames(data);

figure('name','Observables')

for ii=1:numel(vnames)
    
    v=vnames{ii};
    
    subplot(3,2,ii)
    
    plot(data.(v))
    
    title(v)
    
end

%% Normalization

myNormalization={
    'sig_Eps_phi_hetero_1<=sig_Eps_phi_hetero_2'
    };

%% estimate the model
clc

do_estimate=false;

if do_estimate
    
    ms=estimate(m,'data',data,...
        'estim_priors',priors,...
        'kf_presample',3,...
        'kf_init_variance',10,...
        'estim_nonlinear_restrictions',myNormalization);
    
    % save the mode
    %--------------
    
    switch_mode=get(ms,'mode');
    
    save switch_estimates ms switch_mode
    
else
    
    load switch_estimates
    
end

%% irfs

clc

myirfs=irf(ms);

quick_irfs(ms,myirfs)

%% probability plots

plot_probabilities(ms)

%% probability plots

plot_data_against_probabilities(ms)

%% replicate figure 9
close all

probs=ms.filtering.smoothed_state_probabilities;
% I substract a constant since I don't know the correct level of the data.
% Maybe Dan can help with this?
lp=cumsum(data.DLogQl)-0.55;

plotyy(lp,probs.hetero_2,'linewidth',2)

%% (Approximate) Historical decomposition

h=hd(ms);

obslist=get(ms,'obs_list');

figure('name','Observables')

for ii=1:numel(obslist)
    
    v=obslist{ii};
    
    subplot(3,2,ii)
    
    plot_decomp(h.(v))
    
    title(v)
    
    if ii==1
        
        l=legend(h.(v).varnames);
        
    end
    
end

set(l,'interpreter','none')

%% forecast

%% conditional forecast

%% counterfactual

%% posterior sampling
clc

do_sampling=false;

if do_sampling
    
    [objective,lb,ub,mu,SIG]=pull_objective(ms,...
        'solve_check_stability',false,'fix_point_TolFun',1e-6);
    
    sampling_options=struct('MaxTime',3600*24*7,...
        'N',10000,...
        'nchain',4,...
        'thin',5,...
        'adapt_covariance',true);
    
    results=mh_sampler(objective,lb,ub,sampling_options,mu);
    
    save switch_sampling results
    
    % Marginal data density
    algorithms={'mhm','swz','mueller','bridge','is','ris','cj','laplace','laplace_mcmc'};
    options=struct('log_post_kern',objective,'L',10000);
    nalgos=numel(algorithms);
    log_mdd=cell(nalgos,3);
    log_mdd(:,1)=algorithms(:);
    tictoc=nan(nalgos,1);
    
    parfor ii=1:nalgos
        opt=options;
        opt.algorithm=algorithms{ii};
        tic
        log_mdd{ii,2} = mcmc_mdd(results.pop,lb,ub,opt);
        tictoc(ii)=toc;
    end
    log_mdd(:,3)=num2cell(tictoc(:));
    
    save switch_sampling log_mdd -append
    
else
    
    load switch_sampling
    
end

% the benchmark DSGE model (2344.0), the DSGE with two volatility regimes
% (2354.1), and the DSGE model with three volatility regimes (2353.2)
%%

pnames={ms.estimation.priors.name};

drop=0;

obj=mcmc(results.pop,pnames,drop)%,start_from,trimming

%% Trace plots

chain_id=1;

plotfunc=@(name)traceplot(obj,name,chain_id);

hfig=utils.plot.multiple(plotfunc,pnames,'trace plots',3,3)

%% Density plots

chain_id=1;

plotfunc=@(name)densplot(obj,name,chain_id,250);

hfig=utils.plot.multiple(plotfunc,pnames,'density plots',3,3)

%% Mean plots

chain_id=[];

plotfunc=@(name)meanplot(obj,name,chain_id);

hfig=utils.plot.multiple(plotfunc,pnames,'Mean plots',3,3)

%% Autocorrelation plots

chain_id=[];

order=10;

plotfunc=@(name)autocorrplot(obj,name,chain_id,order);

hfig=utils.plot.multiple(plotfunc,pnames,'Autocorrelation plots',3,3)

%% Gelman plots

plotfunc=@(name)gelman_plot(obj,name);

hfig=utils.plot.multiple(plotfunc,pnames,'Gelman (PSRF) plots',3,3)

%% Scatter plots

np=numel(pnames);

scatnames=cell(1,np^2);

iter=0;
for ii=1:np
    for jj=ii+1:np
        iter=iter+1;
        scatnames{iter}={pnames{ii},pnames{jj}};
    end
end
scatnames=scatnames(1:iter);

chain_id=[];

plotfunc=@(name)scatterplot(obj,name{1},name{2},chain_id);%,varargin

hfig=utils.plot.multiple(plotfunc,scatnames,'Scatter plots',7,7)

%% plot priors

[pdata,hdl]=plot_priors(ms,pnames)


%% plot posteriors

plot_posteriors(ms,results.pop,pnames)


%% priors and posteriors

plot_priors_and_posteriors(ms,results.pop,pnames)

%% mode curvature

mode_curvature(ms)
