%  sampler_apt : Shortcut to Adaptive Parallel Tempering algorithm
% 
%  The syntax is ::
% 
%  	results = sampler_apt(target,x0,lb,ub)
% 
%  	results = sampler_apt(target,x0,lb,ub,opts)
% 
%  INPUTS :
% 
%    - *target* : objective to MINIMIZE (Unlike in rsamplers.apt)
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
%        - *H* = 12 : Number of tempering stages
%        - *alpha* = .234 : target acceptance rate
%        - *alpha_swap* = .234 :	target acceptance rate for swaps
%        - *rwm_fixed_p* = 0	: prob of drawing from a fixed initial proposal
%        - *rwm_exp* = 0.6 :	Exponent of random-walk adaptation step size
%        - *sw_exp_rm* = 0.6	: Exp. of temperature adaptation step size
%        - *ram_adapt* = false : Use the robust AM adaptation
%        - *separate_shape_adaptation* = true :separate covariance for each temperature
%        - *fixed_temperatures* = false : fixed temperatures
% 
% 
%  OUTPUT:
% 
%    - *resuls* : cell array of structures
%