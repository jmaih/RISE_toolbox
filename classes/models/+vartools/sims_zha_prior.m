function [abar,SIGu,sampler]=sims_zha_prior(kdata,Yraw,SIGu,...
    prior_hyperparams)

sig=std(Yraw(:,:),[],2); % stretch for panel

ybar=mean(Yraw(:,:),2); % stretch for panel

[Y,X]=vartools.sims_zha_dummies(prior_hyperparams,kdata.nvars,kdata.nx,kdata.nlags,...
    sig,ybar);

% RISE forbids the modification of kdata/self and so we have to create an
% auxiliary

SZkdata=struct('X',kdata.X,'Y',kdata.Y,'linres',kdata.linres);

SZkdata.X=[SZkdata.X(:,:),X];

SZkdata.Y=[SZkdata.Y(:,:),Y];

SZkdata.T=size(SZkdata.Y,2);

[abar,SIGu,sampler]=vartools.sims_zha_posterior(SZkdata.X,...
    SIGu,prior_hyperparams.L5,SZkdata.Y(:),SZkdata.linres);

end