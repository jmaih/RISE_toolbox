function varargout=quantile(this,varargin)
% Computes quantile for time series: It is an interface to MATLAB
%   implementation of quantile 
%
% ::
%
%    varargout = quantile(db,varargin);
%
% Args:
%    db (ts object): times series object
%    varargin: varargin for quantile function in MATLAB
%
% Returns:
%    :
%
%    varargout: output from quantile function
%

% sample quantile (value at %)
[varargout{1:nargout}]=utils.stat.quantile(this,varargin{:});
end