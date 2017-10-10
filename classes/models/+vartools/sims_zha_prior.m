function [abar,SIGu,sampler]=sims_zha_prior(kdata,Yraw,SIGu,...
    prior_hyperparams)

sig=std(Yraw(:,:),[],2); % stretch for panel

ybar=mean(Yraw(:,:),2); % stretch for panel

[Y,X]=vartools.sims_zha_dummies(prior_hyperparams,kdata.nvars,kdata.nx,kdata.nlags,...
    sig,ybar);

kdata.X=[kdata.X(:,:),X];

kdata.Y=[kdata.Y(:,:),Y];

kdata.T=size(kdata.Y,2);

[abar,SIGu,sampler]=vartools.sims_zha_posterior(kdata.X,...
    SIGu,prior_hyperparams.L5,kdata.Y(:),kdata.linres);

end