%  crs_linear_filter_update_cell_estim_K -- filter with update of K. gain
% 
%  ::
% 
% 
%    [loglik,Incr,retcode,Filters]=crs_linear_filter_update_cell_estim_K(...
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
%     - **y** [matrix]: ny x T x npages matrix of data
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
%     - The filter checks for violations of constraints and uses the
%     information in the database to reset the offending variables to their
%     values in the database. When this occurs, the Kalman gain is recomputed
%     so that the traditional updating equation holds.
% 
%     - This strategy is adopted so as to permit the use of the efficient
%     smoothing algorithm of Durbin and Koopman instead of the classical
%     smoothing algorithm which requires multiple inversions of a potentially
%     singular covariance matrix.
% 
%  Example:
% 
%     See also:
%