function varargout=max(varargin)
[varargout{1:nargout}]=utils.stat.nanmax(varargin{:});
end