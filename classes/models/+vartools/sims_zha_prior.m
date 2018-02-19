function [abar,SIGu,sampler]=sims_zha_prior(kdata,Yraw,SIGu,...
    prior_hyperparams)

sig=std(Yraw(:,:),[],2); % stretch for panel

ybar=mean(Yraw(:,:),2); % stretch for panel

[Y,X]=vartools.sims_zha_dummies(prior_hyperparams,kdata.nvars,kdata.nx,kdata.nlags,...
    sig,ybar);

% RISE forbids the modification of kdata/self and so we have to create an
% auxiliary

SZkdata=struct();

SZkdata.estim_.X=[kdata.estim_.X(:,:),X];

SZkdata.estim_.Y=[kdata.estim_.Y(:,:),Y];

SZkdata.estim_.T=size(SZkdata.estim_.Y,2);

SZkdata.estim_.linres=kdata.estim_.linres;

[abar,SIGu,sampler]=vartools.sims_zha_posterior(SZkdata.estim_.X,...
    SIGu,prior_hyperparams.L5,SZkdata.estim_.Y(:),SZkdata.estim_.linres);

end