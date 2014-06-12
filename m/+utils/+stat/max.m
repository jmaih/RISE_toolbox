function varargout=max(varargin)
[varargout{1:nargout}]=nanmax(varargin{:});
end