function varargout=ctw_endogenous_priors(varargin)
% Compute the prior distribution function informed by prior data (past the sampling period)
%
% ::
%
%    [ld, user_info, retcode] = ctw_endogenous_priors(obj);
%    [ld, user_info, retcode] = ctw_endogenous_priors(obj, user_info);
%    [ld, user_info, retcode] = ctw_endogenous_priors(obj, user_info, targets);
%    [ld, user_info, retcode] = ctw_endogenous_priors(obj, user_info, targets, prior_sample);
%
% Args:
%    obj (rise | dsge object):
%    user_info (struct): (default: [])
%    targets (cellstr): names of the target variables
%    prior_sample (int): number of previous sampling period
%
% Returns:
%    :
%
% References:
%    :cite:`christiano2011introducing`
%
% Note:
%    Implemented to match variances
%
% WARNING:
%    As a default, the prior second moments are computed from the same sample
%    used to find the posterior mode. This could be changed by making the
%    appropriate adjustment.
%

[varargout{1:nargout}]=dsge_tools.ctw_endogenous_priors(varargin{:});

end
