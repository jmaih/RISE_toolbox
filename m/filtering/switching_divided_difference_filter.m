%  switching_divided_difference_filter - filter for nonlinear regime-swiching models
% 
%  ::
% 
% 
%    [loglik,Incr,retcode,Filters]=switching_divided_difference_filter(...
%     syst,y,U,z,options)
% 
%  Args:
% 
%     - **syst** [struct]: structure containing:
% 
%           - **PAI00** [vector]: initial probability distributions of regimes
% 
%           - **a** [cell]: initial conditions in each regime
% 
%           - **Qfunc** [function handle]: transition matrix generator
% 
%           - **ff** [function handle]: ft=ff(rt,xt,et), where rt is the
%           regime, xt is the vector of state variables and et the vector of
%           shocks
% 
%           - **P** [cell]: initial covariance matrix of the states in each
%           regime
% 
%           - **H** [cell]: Measurement error covariance matrices in each regime
% 
%           - **SIGeta** [cell]: Covariance matrix of structural shocks.
% 
%     - **y** [matrix]: ny x T matrix of data
% 
%     - **U** [[]|matrix]: ndx x T matrix of exogenous data
% 
%     - **z** [function handle|logical|vector]: linear connection of the
%     observables to the state.
% 
%     - **include_in_likelihood** [logical]: selector of increments to include
%     in the likelihood calculation
% 
%     - **options** [struct]: structure with various options
% 
%  Returns:
%     :
% 
%     - **loglik** [scalar]: log likelihood
% 
%     - **Incr** [vector]: increments of elements going into the likelihood
% 
%     - **retcode** [{0}|integer]: flag for problems.
% 
%     - **Filters** [struct]: Filtered, updated and smoothed variables
% 
%  Note:
% 
%  Example:
% 
%     See also:
%