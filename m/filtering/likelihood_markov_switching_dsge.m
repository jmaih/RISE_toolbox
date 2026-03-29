%  likelihood_markov_switching_dsge : log-likelihood for DSGE models
% 
%  Syntax::
% 
%    [LogLik,Incr,retcode,obj,filtration]=likelihood_markov_switching_dsge(params,obj)
% 
%  Args:
% 
%   - **params** [empty|vector]: parameters under estimation
% 
%   - **obj** [rise|dsge]: model object
% 
%  Returns:
%     
%   - **LogLik** [numeric]: log-likelihood
% 
%   - **Incr** [vector]: density for each period, sum of which gives the
%     log-likelihood
% 
%   - **retcode** [numeric]: return code : 0 if there is no problem
% 
%   - **obj** [rise|dsge]: (possibly) modified model object
% 
%   - **filtration** [struct]: structure with filtering information
% 
% 
%  Note:
% 
%  Example:
% 
%  See also:
%