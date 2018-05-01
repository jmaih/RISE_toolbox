function [is_steady,oldK]=check_steady_state_kalman(is_steady,K,oldK,options,t,no_more_missing)
% check_steady_state_kalman - checks whether the covariance matrix has converged
%
% ::
%
%
%   [is_steady,oldK]=check_steady_state_kalman(is_steady,K,oldK,options,t,no_more_missing)
%
% Args:
%
%    - **is_steady** [true|false]: flag indicating whether or not the steady
%      state is reached
%
%    - **K** [matrix]: Kalman gain
%
%    - **oldK** [matrix]: old Kalman gain
%
%    - **options** [struct]: options for the Kalman filter
%
%    - **t** [integer]: iteration number
%
%    - **no_more_missing** [boolean]: true if there are no more missing
%    observations, false otherwise.
%
% Returns:
%    :
%
%    - **is_steady** [true|false]: flag indicating whether or not the steady
%      state is reached
%
%    - **oldK** [matrix]: updated old Kalman gain
%
% Note:
%
% Example:
%
%    See also:

% the likelihood affects the values of the updated
% probabilities. In turn, the updated probabilities affect the
% values of the predicted probabilities. The predicted
% probabilities enter the collapsing of the covariances and so,
% we never reach the steady state in this case. So far I am
% 100% sure about this. But in order to be 101% sure, I would
% like to run an example.
if no_more_missing && ~is_steady
    discrep=max(abs(K(:)-oldK));
    is_steady=discrep<options.kf_riccati_tol;
    if options.debug
        fprintf(1,'iteration %4.0f  discrepancy %4.4f\n',t,discrep);
    end
end
oldK=K(:);
end
