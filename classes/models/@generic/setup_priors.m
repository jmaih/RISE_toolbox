%--- help for generic/setup_priors ---
%
%  INTERNAL FUNCTION: Format the priors for the estimation of models in RISE.
% 
%  ::
% 
%    obj=setup_priors(obj,MyPriors)
%    obj=setup_priors(obj,MyPriors,error_control)
% 
%  Args:
% 
%     obj (rise | dsge | svar | rfvar): model object. This function is not
%       meant to be called directly by the user. But it shows how to write the
%       priors when those are not written inside the model file.
% 
%     MyPriors (struct): Each field of the structure is one of the
%       following alternatives
% 
%       - Maximum likelihood or uniform priors::
% 
%           P.pname={start_value,lower_bound,upper_bound};
% 
%       - Bayesian prior using mean and standard deviation::
% 
%           P.pname={start_value,prior_mean,prior_stdev,'distribution'};
% 
%       - Same as above except for adding a hard lower bound::
% 
%           P.pname={start_value,prior_mean,prior_stdev,'distribution',lower_bound};
% 
%       - Same as above except for adding a hard upper bound::
% 
%           P.pname={start_value,prior_mean,prior_stdev,'distribution',lower_bound,upper_bound};
% 
%       - Bayesian prior using quantiles of the distribution::
% 
%           P.pname={start_value,lower_quantile,lower_quantile,'distribution(prob)'};
% 
%       - Same as above except for adding a hard lower bound::
% 
%           P.pname={start_value,lower_quantile,lower_quantile,'distribution(prob)',lower_bound};
% 
%       - Same as above except for adding a hard upper bound::
% 
%           P.pname={start_value,lower_quantile,lower_quantile,'distribution(prob)',lower_bound,upper_bound};
% 
%       - Dirichlet
% 
%           Denote by dirichlet_j, the jth dirichlet distribution, j=1,2,... we
%           have P.dirichlet_j={sd_ii,pname1,m1,pname2,m2,...,pnamen,mn} where
% 
%             - sd_ii: is the standard deviation of the diagonal element (a
%               parameter never listed by RISE) of the dirichlet distribution.
%             - pname1, pname2,...,pnamen: are the names of the off diagonal
%               parameters
%             - m1, m2,...,mn: are the means of each parameters
% 
%     error_control (empty | cell): element constructed by RISE for
%       controling the syntax used in the model file.
% 
%  Returns:
%     :
% 
%     - **obj** [rise|dsge|svar|rfvar]: updated model object.
% 
%  Note:
% 
%     - This function is also indirectly used for svar and rfvar objects.
%     - Possible distributions include:
% 
%        - beta
%        - cauchy
%        - gamma
%        - inv_gamma
%        - laplace
%        - left_triang
%        - logistic
%        - lognormal
%        - normal
%        - pareto
%        - right_triang
%        - uniform
%        - weibull
% 
%  See also:
%     - rise_generic/estimate
%     - dsge/estimate
% 
%