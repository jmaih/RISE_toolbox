%--- imh.m not found. Showing help for imm instead. ---
%
%  IMM performs the constrained regime-switching Kalman filter for state
%  space models using the interactive Multiple Model algorithm of Blom and
%  Bar-Shalom (1988)
% 
%    [loglik] = imm(syst, data_y, U, z, options)
%    [loglik, Incr, retcode, Filters] = imm(...)
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
%        [loglik, Incr, retcode, Filters] = imm(syst, data_y, U, z, options);
% 
%    See also:
% 
%    Reference: Blom and Bar-Shalom (1998): "The Interacting Multiple Model
%    Algorithm for Systems with Markovian Switching Coefficients"
%