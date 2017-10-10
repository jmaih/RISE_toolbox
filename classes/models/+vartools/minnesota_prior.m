function [abar,SIGu,sampler]=minnesota_prior(kdata,Yraw,SIGu,...
    prior_hyperparams)

astar=vartools.set_prior_mean(kdata,prior_hyperparams);

Va = vartools.set_prior_variance(Yraw,SIGu,kdata,prior_hyperparams);

[abar,SIGu,sampler]=vartools.normal_prior_with_known_Sigma(kdata.X,...
    SIGu,Va,astar,kdata.Y(:),kdata.linres);

end