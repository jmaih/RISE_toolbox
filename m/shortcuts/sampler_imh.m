%  sampler_imh : Shortcut to Independent Metropolis-Hastings algorithm
% 
%  The syntax is ::
% 
%  	results = sampler_imh(target,x0,lb,ub)
% 
%  	results = sampler_imh(target,x0,lb,ub,opts)
% 
%  INPUTS :
% 
%    - *target* : objective to MINIMIZE (Unlike in rsamplers.imh)
%    - *x0* : vector of initial conditions
%    - *lb* : lower bound
%    - *ub* : upper bound
%    - *opts* : options. These include:
% 
%        - *c* = 1 : cov scaling parameter
%        - *tunedCov* : covariance matrix of the parameters
%        - *proposal* = 'normal' :	proposal in {'normal','t-student'}
%        - *nchain* = 1 : number of chains
%        - *N* = 2000 : number of draws
%        - *thinning* =1 : number of thinning draws
%        - *burnin* = 0 : burnin sample
%        - *MaxTime* = inf : Max Time
%        - *MaxFunEvals* = inf :
% 
% 
%  OUTPUT:
% 
%    - *resuls* : cell array of structures
%