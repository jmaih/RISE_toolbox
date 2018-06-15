function varargout=quantile(this,varargin)
% INTERNAL FUNCTION
%

% sample quantile (value at %)
[varargout{1:nargout}]=utils.stat.quantile(this,varargin{:});
end