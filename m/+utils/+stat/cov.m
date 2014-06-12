function varargout=cov(varargin)
[varargout{1:nargout}]=nancov(varargin{:});
end