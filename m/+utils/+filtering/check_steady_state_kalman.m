%  INTERNAL FUNCTION: checks whether the covariance matrix has converged
% 
%  ::
% 
%    [is_steady,oldK]=check_steady_state_kalman(is_steady,K,oldK,options,t,no_more_missing)
% 
%  Args:
% 
%     - **is_steady** [true|false]: flag indicating whether or not the steady
%       state is reached
%     - **K** [matrix]: Kalman gain
%     - **oldK** [matrix]: old Kalman gain
%     - **options** [struct]: options for the Kalman filter
%     - **t** [integer]: iteration number
%     - **no_more_missing** [boolean]: true if there are no more missing
%       observations, false otherwise.
% 
%  Returns:
%     :
% 
%     - **is_steady** [true|false]: flag indicating whether or not the steady
%       state is reached
%     - **oldK** [matrix]: updated old Kalman gain
% 
%