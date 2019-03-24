%--- help for generic/log_posterior_kernel ---
%
%  Computes the log posterior of the dsge model
% 
%  Note:
%     In effort to make RISE modular, this function is available so that one
%     can use a different sampler if needed, but most likely, one should
%     just use available stock samplers.  
% 
%  ::
% 
%     [log_post,log_lik,log_prior,Incr,retcode,obj]=log_posterior_kernel(obj, param)
% 
%  Args:
%     model (dsge | rise object): model object
%     param (column vector): parameter values
% 
%  Returns:
%     : [log_post, log_lik, log_prior, Incr, retcode, model]
% 
%     - **log_post** (double): log posterior
%     - **log_lik** (double): log likelihood
%     - **log_prior** (double): log prior
%     - **Incr** (double):
%     - **retcode** : return code
%     - model (dsge | rise object): model passed through
% 
%