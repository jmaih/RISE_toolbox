%--- help for abstvar/autocov ---
%
%  Compute the autocovariances (and the auto-correlation) of endogenous variables given the parameter values
% 
%  ::
% 
%     [C,R] = autocov(self);
%     [C,R] = autocov(self, params);
%     [C,R] = autocov(self, params, max_periods);
% 
%  Args:
%     self (var object): var object
%     params (cell of struct): struct containing var model related parameters (default: [])
%     max_periods (integer): maximum number of period to calculate auto-covariance (default: 5)
% 
%  Returns:
%     :
% 
%     - **C** [4-dimensional array]: auto-covariance where dimensions correspond to
% 
%        - 1,2: covariance
%        - 3: time lags
%        - 4: Number of parameters
% 
%     - **R** [4-dimensional array]: auto-correlation with same dimensions
% 
%