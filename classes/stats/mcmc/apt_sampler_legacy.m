%  MH_SAMPLER -- Metropolis Hastings sampler
% 
%  ::
% 
% 
%    [Results]=MH_SAMPLER(logf,lb,ub)
% 
%    [Results]=MH_SAMPLER(logf,lb,ub,options)
% 
%    [Results]=MH_SAMPLER(logf,lb,ub,options,mu)
% 
%    [Results]=MH_SAMPLER(logf,lb,ub,options,mu,SIG)
% 
%  Args:
% 
%     - **logf** [char|function_handle]: Objective function to MINIMIZE!!!
% 
%     - **lb** [d x 1 vector]: lower bound of the paramters
% 
%     - **ub** [d x 1 vector]: upper bound of the paramters
% 
%     - **options** [struct]:
%       - **alpha** [scalar|2-element|{[.25,.45]}]: target acceptance rate
%       - **burnin** [integer|{0}]: number of burn-in initial simulations
%       - **N** [integer|{20000}]: number of simulations
%       - **verbose** [integer|{100}]: displays progress for every multiple of
%       "verbose"
%       - **c** [scalar|{0.25}]: initial scale for the covariance matrix
%       - **c_range** [vector|{}]: range of variation of c
%       - **thin** [integer|{1}]: number of thinning simulations. 1 means we
%       keep every draw, 2 means we keep every second, 3 every third, etc.
%       - **retune_cov_every** [integer|{100}]: frequence for the retuning of
%       the scale parameter
%       - **penalty** [numeric|{1e+8}]: worst possible function value
%       - **nchain** [integer|{1}]: number of parallel chains
%       - **rwm_exp** [numeric|{0.6}]: tuning hyper-parameter for scale and
%       covariance matrix
%       - **fixed_scaling** [true|{false}]: if true, the scaling (c) of the
%       covariance matrix is kept constant
%       - **use_true_moments** [true|{false}]: if true, the updated exact
%       covariance matrix of the draws is used at each step. If false, a
%       different update of the covariance matrix is used.
%       - **logproppdf** [function_handle|{[]}]: used when the proposal is not
%       symmetric
%       - **MaxTime** [numeric|{inf}]: maximum simulation time.
%       - **adapt_covariance** [true|{false}]: If true, the covariance matrix
%       is updated with the sampled draws.
%       - **save** [struct|{'every='inf,'location=',pwd,'filename=',''}]:
%       structure with fields "every", "location", "filename"
%       - **recover** [false|{true}]: attempts to recover from an earlier
%       aborted simulation
%       - **recover_start_at_best** [true|{false}]: in an event of a recovery,
%       start all the chains at the best value. Note this will break the chains
%       if the best value is not the last element that was saved in the chain.
% 
%     - **mu** [d x 1 vector]: initial condition for the sampler
% 
%     - **SIG** [d x d matrix]: initial covariance matrix
% 
%  Returns:
%     :
% 
%     - **Results** [struct]:
%       - **pop** [nchain x N struct]: fields are "x" for the parameter vector
%       and "f" for the value of the parameter vector
%       - **bestf** [numeric]: best function value
%       - **bestx** [vector]: best parameter vector
%       - **best** [nchain x 1]: vector of best individual in each chain
%       - **m** [vector]: mean of the parameter draws
%       - **SIG** [matrix]: covariance of the parameter draws
%       - **m_algo** [vector]: mean with particular updating
%       - **SIG_algo** [matrix]: covariance with particular updating
%       - **funevals** [1 x nchain vector]: function evaluations
%       - **stats** [struct]: stats on the optimization
% 
%  Note:
% 
%     - It is assumed that logf is a function to minimize
% 
%  Example:
% 
%     See also:	CONSTANT_BVAR_SAMPLER, RRF_SAMPLER
%