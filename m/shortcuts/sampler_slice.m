%  sampler_slice : Shortcut to Slice algorithm
% 
%  The syntax is ::
% 
%  	results = sampler_slice(target,x0,lb,ub)
% 
%  	results = sampler_slice(target,x0,lb,ub,opts)
% 
%  INPUTS :
% 
%    - *target* : objective to MINIMIZE (Unlike in rsamplers.slice)
%    - *x0* : vector of initial conditions
%    - *lb* : lower bound
%    - *ub* : upper bound
%    - *opts* : options. These include:
% 
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