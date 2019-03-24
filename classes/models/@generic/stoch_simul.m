%--- help for generic/stoch_simul ---
%
%  stoch_simul -- attempts to emulate dynare's stoch_simul
% 
%  ::
% 
%     oo_=stoch_simul(obj)
% 
%     oo_=stoch_simul(obj,var_list)
% 
%     oo_=stoch_simul(obj,var_list,varargin)
% 
%  Args:
% 
%     obj (dsge | rise ): model object
% 
%     var_list (empty | char | cellstr ): List of variables for which to run
%       stoch_simul
% 
%     varargin : Pairwise list of extra arguments
% 
%  Returns:
%     :
% 
%     - **oo_** [struct]: structure containing
%        - **irfs** [struct]: structure containing impulse responses in time
%          series format
%        - **simulations** [struct]: structure containing impulse responses 
%          in time series format
%        - **vcov** [matrix]: (simulated) variance-covariance
%        - **skewness** [vector]: skewness measure
%        - **kurtosis** [vector]: kurtosis measure
%        - **variance** [vector]: (simulated) variances for the endogenous
%          variables 
%        - **stdev** [vector]: (simulated) standard deviations for the endogenous
%          variables 
%        - **corrcoef** [matrix]: (simulated) correlation coefficients for
%          the endogenous variables
% 
%