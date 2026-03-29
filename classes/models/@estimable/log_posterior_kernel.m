%--- help for estimable/log_posterior_kernel ---
%
%  Computes the log posterior of the dsge model
% 
%  ::
% 
%     [log_post,log_lik,log_prior,Incr,retcode,obj]=log_posterior_kernel(obj, param)
% 
%  Args:
% 
%     - obj (estimable object): model object
% 
%     - param (column vector): parameter values
% 
%  Returns:
% 
%     - **log_post** (double): log posterior
% 
%     - **log_lik** (double): log likelihood
% 
%     - **log_prior** (double): log prior
% 
%     - **Incr** (double):
% 
%     - **retcode** : return code
% 
%     - obj (estimable object): model passed through
% 
%  See also:
% 
%     - :ref:`log_prior_density <rise.log_prior_density>`
% 
% 
%  .. note:: In effort to make RISE modular, this function is available so that one
%     can use a different sampler if needed, but most likely, one should
%     just use available stock samplers.
%