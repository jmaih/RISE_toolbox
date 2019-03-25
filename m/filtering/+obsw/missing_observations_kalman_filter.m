%  INTERNAL FUNCTION: Kalman filter for non-regime-switching models
% 
%  ::
% 
%    f=missing_observations_kalman_filter(y,T,R,Z,H,Q,init,dy,ca,level)
% 
%  Args:
% 
%       y (p x n matrix): data
%       T (m x m x nT): state matrix, see below
%       R (m x r x nR): state matrix, see below
%       Z (p x m | logical): state matrix or selector, see below
%       H (p x p x nH): Covariances of measurement errors
%       Q (r x r x nQ): Covariances of shocks
%       init (struct):
%       dy (p x ndy):
%       ca (m x nca):
%       level (0 | 1 | 2 | {3}): controls the level output produced.
% 
%           - if level==0, only the likelihood is returned
%           - if level==1, filters/forecasts are returned in addition to level 0
%           - if level==2, updates are returned in addition to level 1
%           - if level==3, smoothed values are returned in addition to level 2
% 
%       nstep (integer | {1}): Number of forecast steps
% 
%  Returns:
%     :
% 
%       - **f** [struct]: output structure
% 
%           - **retcode** [0|305]: 305 if computation of the likelihood breaks
%             downs
%           - **incr** [vector]: likelihood in each period
%           - **log_lik** [scalar]: log likelihood
%           - **a** [m x n+1 x nstep]: filtered stated (available if level >0)
%           - **P** [m x m x n+1 x nstep]: Covariance matrix of filtered state
%             (available if level >0)
%           - **att** [m x n]: Updated state vector  (available if level >1)
%           - **Ptt** [m x m x n]: Covariance matrix for updates (available if
%             level >1)
%           - **iF** [p x p x n]: Inverses of covariance matrix of forecast
%             errors
%           - **v** [p x n]: Forecast errors
%           - **K** [m x p x n]: Kalman gains
%           - **alpha** [m x n]: smoothed state vector
%           - **V** [m x m x n]: Covariance matrix of smoothed state vector
%           - **r** [m x n+1]: quantity useful for the efficient computation of
%             the state vector
%           - **eta** [r x n]: smoothed shocks
%           - **epsilon** [p x n]: smoothed measurement errors
% 
%  Note:
% 
%     The system is of the form
%     .. math::
% 
%        alpha(t+1) = T*alpha(t) + ca(t) + R(t)*eta(t),  eta ~ N(0, Q)
%        y(t) =   Z*alpha(t) + dy(t) + epsilon(t),  epsilon ~ N(0,H)
% 
%