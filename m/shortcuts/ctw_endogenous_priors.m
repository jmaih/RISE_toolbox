function varargout=ctw_endogenous_priors(varargin)

% Endogenous priors a la 
% Christiano,Trabandt and Walentin (2011):"Introducing Financial Frictions
% and Unemployment into a Small Open Economy Model", Journal of Economic
% Dynamics and Control
% the priors include a metric across some choosen moments of the (supposedly
% pre-sample) data.
% *** Implemented file for variances, but in principle any moment
% *** could be matched
% As a default, the prior second moments are computed from the same sample
% used to find the posterior mode. This could be changed by making the
% appropriate adjustment.

[varargout{1:nargout}]=dsge_tools.ctw_endogenous_priors(varargin{:});

end
