function [abar,SIGu,sampler]=normal_wishart_prior(kdata,Yraw,SIGu,...
    prior_hyperparams)

astar=vartools.set_prior_mean(kdata,prior_hyperparams);

[~,V] = vartools.set_prior_variance(Yraw,SIGu,kdata,prior_hyperparams);

[abar,SIGu,sampler]=vartools.normal_wishart_posterior(kdata.estim_.X,...
    SIGu,V,astar,kdata.estim_.Y(:),kdata.linres);

end