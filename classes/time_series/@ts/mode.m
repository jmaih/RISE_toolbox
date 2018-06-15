function varargout=mode(this,varargin)
% INTERNAL FUNCTION
%

nout=nargout;
[varargout{1:nout}]=utils.stat.mode(this.data,varargin{:});
end
