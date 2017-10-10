function bstar=set_prior_mean(kdata,prior_hyperparams)

coefprior=prior_hyperparams.coefprior;

B=vartools.bvar_coef_prior(coefprior,kdata.nvars,kdata.nlags,kdata.nx);

bstar=B(:);

end