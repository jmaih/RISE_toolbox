%--- help for abstvar/irf ---
%
%  Compute impulse response function from the given parameter values
% 
%  ::
% 
%     myirfs = irf(self, shock_names, irf_periods, params, Rfunc);
% 
%  Args:
% 
%     self (var object): var object
% 
%     shock_names (cellstr): shocks to compute IRFs (has to be consistent with the names given in :func:`identification <var.identification>`)
% 
%     irf_periods (int): number of periods to compute IRFs (default: 40)
% 
%     params (vector): parameter values of the model. If empty MLE/posterior-mode values are used.
% 
%     Rfunc (function handle): identification function. This is an output of :func:`identification <var.identification>`. (default: choleski identification)
% 
%     girf_setup (struct|{empty}): structure containing information relevant
%        for the computation of generalized impulse response functions. If
%        empty, simple regime-specific impulse responses are computed, else
%        girfs are computed. In that case the relevant information to
%        provide in girf_setup is:
%        - nsims : (default=300) number of simulations for the integration.
%        Note that even setting girf_setup=struct() will trigger the
%        computation of girfs. But in that case only the default options
%        will apply.
% 
%  Returns:
% 
%     : struct containing IRFs
% 
%  Note:
% 
%     Only successful IRFs are returned. If the structure does not return some IRFs make sure that IRFs properly identified.
% 
%
%    Other functions named irf
%
%       dsge/irf    generic/irf
%