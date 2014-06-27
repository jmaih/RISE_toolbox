%% Test of distributions
clc
close all
clear all
%% univariate distributions
%                name, start, plb, pub
Distributions={
    'beta',0.5,0.05,0.948
    'normal',0.5,-2,2 % -inf < x < inf
    'inv_gamma',0.5,0.005,1% 0 < x < inf
    'lognormal',1.5,0.5,3 % 0 < x < inf
    'truncated_normal',0.5,0.05,0.948
    'uniform',0.5,-1,3 % -inf < x < inf
    'weibull',1.5,0.5,3 % 0 < x < inf
    'pareto',3.5,3,7  % 0 < x < inf
    'cauchy',0.5,-20,30  % -inf < x < inf
    'laplace',0.5,-1,3 % -inf < x < inf
    'logistic',0.5,-1,3 % -inf < x < inf
    'gamma',1.5,0.5,7  % 0 < x < inf
    };
                        
obj=rise_estim_param.empty(0);
for ii=1:size(Distributions,1)
    distrib=Distributions{ii,1};
    name=['par_',distrib];
    tex_name='junk';
    value=Distributions{ii,2};
    plb=Distributions{ii,3};
    pub=Distributions{ii,4};
    prob=0.9;
    prior_trunc=1e-10;
    obj(ii,1)=rise_estim_param(name,tex_name,ii,value,plb,pub,distrib,prob,prior_trunc);
end
%% plot it all
plot(obj)
%% multivariate: mv_normal, inv_wishart, wishart, dirichlet, truncated_mv_normal
% support: hyperparameter_residuals find_hyperparameters  
