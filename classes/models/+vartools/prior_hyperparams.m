function p=prior_hyperparams()

% L1 : Overall tightness
% L2 : Cross-variable specific variance parameter
% L3 : Speed at which lags greater than 1 converge to zero
% L4 : tightness on deterministic/exogenous terms
% L5 : covariance dummies(omega)
% L6 : co-persistence (lambda)
% L7 : Own-persistence (mu)

% prior types
% 'minnesota' (default)
% 'normal-wishart','nw'
% 'indep-normal-wishart','inw'
% 'jeffrey'
% 'sims-zha','sz'

p=struct();

p.L1=0.1; %0.5;

p.L2=0.5;

p.L3=1; % 2

p.L4=10^2;

p.coefprior=[];

p.type='minnesota';

p.normal_wishart_eta=0.01;

% a small value of ? implies a tight prior variance
% If ? is increased, the posterior mean approaches the LS estimator as
% expected because all terms involving V ?1 in the formulas for the
% posterior moments disappear if ??? and V ?1 ? 0 

p.independent_normal_wishart_eta=p.normal_wishart_eta;

% further hyperparameters for dummy observations

% covariance dummies(omega)
p.L5=3;

% co-persistence (lambda)
p.L6=5;

% Own-persistence (mu)
p.L7=2;

end