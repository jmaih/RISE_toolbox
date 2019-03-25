%  Initial draws for populations mcmc algorithms
% 
%  ::
% 
%    [pop,funevals,log_I_init,covar]=initial_draws(logf,lb,ub,N)
%    [pop,funevals,log_I_init,covar]=initial_draws(logf,lb,ub,N,penalty)
%    [pop,funevals,log_I_init,covar]=initial_draws(logf,lb,ub,N,penalty,mu)
% 
%  Args:
% 
%     - **logf** [function handle]: objective to **minimize**
%     - **lb** [vector]: lower bound
%     - **ub** [vector]: upper bound
%     - **N** [integer]: number of parameters to draws
%     - **penalty** [numeric|{1e+8}]: maximum acceptable value of the objective function
%     - **mu** [vector]: start draws
%     - **max_attempts** [integer|{50}]: maximum number of attempts for each
%       parameter vector.
% 
%  Returns:
%     :
% 
%     - **pop** [struct]: initial draws, each with fields
% 
%       - **f** [numeric]: value of fitness
%       - **x** [vector]: parameters
% 
%     - **funevals** [integer]: number of function evaluations
%     - **log_I_init** [0]: initial value of the log marginal data density for
%       a multivariate uniform distribution
%     - **covar** [matrix]: variance-covariance matrix of the drawn parameters
%     - **NumWorkers** [integer]: number of workers available for parallel
%       processing.
% 
%