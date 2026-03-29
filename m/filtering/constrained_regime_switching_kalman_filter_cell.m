% CONSTRAINED_REGIME_SWITCHING_KALMAN_FILTER_CELL performs the constrained
% regime-switching Kalman filter for state space models. 
% 
%    [loglik] = constrained_regime_switching_kalman_filter_cell(syst, data_y, U, z, options)
%    [loglik, Incr, retcode, Filters] = constrained_regime_switching_kalman_filter_cell(...)
% 
%    Args:
%        - `syst` : State space model structure.
%        - `data_y` : Observations.
%        - `U` : Exogenous inputs.
%        - `z` : Measurement equation structure.
%        - `options` : Options structure.
% 
%    Returns:
%        - `loglik` : Log-likelihood of the model.
%        - `Incr` : Incremental likelihood.
%        - `retcode` : Return code indicating success (0) or an error.
%        - `Filters` : Structure containing filtering results.
% 
%    This function performs the constrained regime-switching Kalman filter
%    for state space models with a specific structure. 
% 
%    Example:
%        [loglik, Incr, retcode, Filters] = constrained_regime_switching_kalman_filter_cell(syst, data_y, U, z, options);
% 
%    See also: 
% 
%    Reference:
%