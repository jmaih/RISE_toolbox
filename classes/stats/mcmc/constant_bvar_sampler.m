%  CONSTANT_BVAR_SAMPLER -- Sampler for constant-parameter BVARs with or
%  without parameter restrictions
% 
%  ::
% 
% 
%    [Results]=CONSTANT_BVAR_SAMPLER(start)
% 
%    [Results]=CONSTANT_BVAR_SAMPLER(start,options)
% 
%  Args:
% 
%     - **start** [struct]: structure containing the fields of interest for
%     posterior sampling:
%       - **prior_type** ['diffuse'|'jeffrey'|'minnesota'|'normal_wishart'|...
%           'indep_normal_wishart']: type of prior
%       - **Y** [matrix]: Left-hand side
%       - **X** [matrix]: right-hand side
%       - **nobs** [numeric]: number of observations
%       - **K** [numeric]: number of parameters per equation
%       - **a_func** [function_handle]: function that inflates x and Vx under
%           linear restrictions.
%       - **a2tilde** [struct]:
%           - **prior** []:
%           - **ols** []:
%           - **post** []:
%       - **na2** [numeric]: number of unrestricted parameters
%       - **estimafy** [function_handle]: computes posterior and ols
%       - **a2tilde_func** [function_handle]: shortens the vector of parameters
%       to estimate
% 
%     - **options** [struct]:
%       - **burnin** [integer|{0}]: number of burn-in initial simulations
%       - **N** [integer|{20000}]: number of simulations
%       - **thin** [integer|{1}]: number of thinning simulations. 1 means we
%       keep every draw, 2 means we keep every second, 3 every third, etc.
%       - **MaxTime** [numeric|{inf}]: maximum simulation time.
% 
%  Returns:
%     :
% 
%     - **Results** [struct]:
%       - **pop** [1 x N struct]: fields are "x" for the parameter vector
%       and "f" for the value of the parameter vector
%       - **m** [vector]: mean of the parameter draws
%       - **SIG** [matrix]: covariance of the parameter draws
% 
%     - **start** [struct]: see above
% 
%  Note:
% 
%  Example:
% 
%     See also: MH_SAMPLER
%