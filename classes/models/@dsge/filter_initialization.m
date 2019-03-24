%--- help for filter_initialization ---
%
%  Initial conditions for filtering
% 
%  ::
% 
%    [init,retcode]=filter_initialization(obj)
% 
%  Args:
% 
%     obj (rise | dsge): model object
% 
%     varargin (name,value): valid pairwise options with the most
%       relevant beeing:
% 
%       - **kf_ergodic** [{true}|false]: initialization at the ergodic
%         distribution
%       - **kf_init_variance** [{[]}|scalar]: initial variance factor (Harvey
%         scale factor). If not empty, the information in T and R is ignored.
%       - **kf_presample** [{0}|integer]: Number of observations to discard
%         before computing the likelihood.
%       - **kf_filtering_level** [0|1|2|{3}]: 0: Likelihood only, 1: 0+filtered
%         series, 2: 1+ updated series, 3: 2+ smoothed series
%       - **kf_user_init** [{[]}|cell]: User-defined initialization. When not
%         empty, it can take three forms. {a0}, {a0,cov_a0}, {a0,cov_a0,PAI00}
%         where a0 is the initial state vector with the same order as the rows of
%         T, cov_a0 is the initial covariance of the state vector (same order as
%         a0) and PAI00 is the initial vector of regime probabilities.
%       - **kf_user_algo** [{''}|char|function handle]: User-defined filtering
%         algorithm. It should have the same inputs and outputs as e.g.
%         switching_divided_difference_filter.m.
%       - **kf_householder_chol** [{false}|true]: if true, return the cholesky
%         decomposition when taking the householder transformation. This option
%         is primarily used in the switching divided difference filter.
% 
%  Returns:
%     :
% 
%     - **init** [struct]: initial conditions of the filter
%     - **retcode** [scalar]: 0 if there is no problem
%     - **nsols** [scalar]: number of solutions
% 
%  Note:
% 
%     - The initial covariance is always computed using the current impact
%       matrix of the shocks even when the anticipation horizon of the agents is
%       greater than zero.
%     - use option "lyapunov_diffuse_factor" to set the variance for
%       nonstationary variables separately. If this is not the case the variance
%       for those variables will be 1 by default.
% 
%  See also:
%     - dsge/filter
% 
%