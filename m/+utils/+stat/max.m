function varargout=max(varargin)
% INTERNAL FUNCTION
%

[varargout{1:nargout}]=utils.stat.nanmax(varargin{:});
end