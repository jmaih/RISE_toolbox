%--- help for generic/log_prior_density ---
%
%  Computes the probability density function of the prior corresponding to
%  the parameter values
% 
%  Note:
%     In effort to make RISE modular, this function is available so that one
%     can use a different sampler if needed, but most likely, one should
%     just use available stock samplers.  
% 
%  ::
% 
%     [lnprior, retcode] = log_prior_density(model)
%     [lnprior, retcode] = log_prior_density(model, param)
%     [lnprior, retcode] = log_prior_density(model, param,filtration)
% 
%  Args:
%     model (dsge | rise object): model object
% 
%     param (column vector): parameter values
% 
%     filtration (empty|struct): results from model filtration that can
%     potentially be used for endogenizing priors.
% 
% 
%  Returns:
%     : [lnprior, retcode]
% 
%     - **lnprior** (double): log of prior density function
%     - **retcode**: return code
% 
%  See also:
%     - :func:`log_posterior_kernel <dsge.log_posterior_kernel>`
% 
%