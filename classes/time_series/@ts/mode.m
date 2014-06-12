function varargout=mode(this,varargin)
nout=nargout;
[varargout{1:nout}]=utils.stat.mode(this.data,varargin{:});
end
