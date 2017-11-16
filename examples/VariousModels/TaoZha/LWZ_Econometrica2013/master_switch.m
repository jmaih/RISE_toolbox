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

ms=estimate(m,'data',data,...
    'estim_priors',priors,...
    'kf_presample',3,...
    'kf_init_variance',10,...
    'estim_nonlinear_restrictions',myNormalization);

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

[objective,lb,ub,mu,SIG]=pull_objective(ms,...
    'solve_check_stability',false,'fix_point_TolFun',1e-6);

sampling_options=struct('MaxTime',3600*24*7,...
    'N',10000,...
    'nchain',4,...
    'thin',5,...
    'adapt_covariance',true);

results=mh_sampler(objective,lb,ub,sampling_options,mu);

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
%% priors and posteriors

%% mode curvature


