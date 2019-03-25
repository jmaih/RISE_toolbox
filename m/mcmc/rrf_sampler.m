%  RRF_SAMPLER -- Rapid Reaction Force Sampler
% 
%  ::
% 
% 
%    [Results]=RRF_SAMPLER(logf,lb,ub)
% 
%    [Results]=RRF_SAMPLER(logf,lb,ub,options)
% 
%    [Results]=RRF_SAMPLER(logf,lb,ub,options,mu)
% 
%    [Results]=RRF_SAMPLER(logf,lb,ub,options,mu,SIG)
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
%       - **c** [scalar|{1}]: initial scale for the covariance matrix
%       - **c_range** [vector|{[sqrt(eps),100]}]: range of variation of c
%       - **alpha** [scalar|2-element|{[.25,.45]}]: target acceptance rate
%       - **burnin** [integer|{0}]: number of burn-in initial simulations
%       - **N** [integer|{20000}]: number of simulations
%       - **M_over_N** [numeric|{5/100}]: percentage of draws per bin
%       - **H** [integer|{40}]: number of tempering stages
%       - **extra_runs** [integer|{5}]: number of extra runs after lambda=1
%       - **lambda_1** [numeric|{2.5e-5}]: initial tempering
%       - **ps** [numeric|{0.9}]: probability of storing a draw
%       - **p** [numeric|{[]}]: probability of drawing from previous
%       distribution
%       - **p_mutant** [numeric|{0}]: probability of drawing a mutant
%       - **ess_min** [numeric|{0.1}]: constant for controlling the evolution
%       of lambda
%       - **focus_fire_after** [numeric|{0.5}]: percentage after which
%       exploitation takes over exploration
%       - **penalty** [numeric|{1e+8}]: worst possible function value
%       - **fixed_scaling** [true|{false}]: if true, the scaling (c) of the
%       covariance matrix is kept constant
%       - **geometric_lambda** [{true}|false]: geometric evolution of lambda
%       - **rwm_exp** [numeric|{0.6}]: tuning hyper-parameter for scale and
%       covariance matrix
%       - **use_true_moments** [true|{false}]: if true, the updated exact
%       covariance matrix of the draws is used at each step. If false, a
%       different update of the covariance matrix is used.
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
%       - **lambda** [vector]: tempering levels
%       - **best** [nchain x 1]: vector of best individual in each chain
%       - **m** [vector]: mean of the parameter draws
%       - **SIG** [matrix]: covariance of the parameter draws
%       - **m_algo** [vector]: mean with particular updating
%       - **SIG_algo** [matrix]: covariance with particular updating
%       - **funevals** [integer]: function evaluations
%       - **log_mdd** [numeric]: log marginal data density
%       - **ess** [numeric]: effective sample size
%       - **individual_cores_stats** [struct]: stats on the optimization
% 
%  Note:
% 
%     - It is assumed that logf is a function to minimize
% 
%     - the update of c in one bin might not be adequate for the other bins. As
%     a result, the last bin may be detrimental for the first bin. This is
%     calls for bin-specific log_c: IMPLEMENTED!!!
% 
%     - Possible variations for the stud
%       - best in the bin (Current implementation)
%       - random
%       - simple average : what if it has unacceptable density?
%       - weigthed average through a merit system : what if it has unacceptable density?
% 
%  Example:
% 
%     See also:	CONSTANT_BVAR_SAMPLER, MH_SAMPLER
%