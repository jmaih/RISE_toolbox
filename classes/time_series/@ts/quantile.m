function varargout=quantile(this,varargin)

% sample quantile (value at %)
[varargout{1:nargout}]=utils.stat.quantile(this,varargin{:});
end