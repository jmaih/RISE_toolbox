%% housekeeping
close all
clear
clc
%% load the data

load tut01_data

%% set up VAR model

warning('This system is nonstationary because of LGDP !!!')

endog={'rtwi','LGDP','Aldp','r','LRER'};

nlags=4;

const=true;

exog={'du93Q1','du95Q4','du92Q3'};

v=rfvar(endog,exog,nlags,const);

%%  set prior

var_prior=rfvar.prior_template();
% modify as needed

var_prior.type='sz';

prior=struct('var',var_prior);%,'nonvar',switch_prior

%%
ve=estimate(v,db,{'1982Q1','2004Q4'},prior);%,restrictions

%% should the Feds fund rate react to domestic variables?
linres={};

for ilag=1:nlags
    
    for iv=2:numel(endog)
        
        y=endog{iv};
        
        linres=[linres;{sprintf('b%0.0f(1,%s)=0',ilag,y)}];
        
    end
    
end

%%
ve_lr=estimate(v,db,{'1982Q1','2004Q4'},prior,linres);

%% sampler

% tmp=ve.estim_.sampler(200)
% tmp=ve_lr.estim_.sampler(200)
%
% N.B: samples returned in terms of untransformed restrictions i.e. original
% parameters

%% save estimated models for later use

models=struct('ve',ve,'ve_lr',ve_lr);

save('tut13_bayestimation','models')
									
