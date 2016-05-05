function varargout=dynare2rise(varargin)
% converts a basic dynare file into a RISE one
%
% dynare2rise(dynFileName,riseFileName)
% dynare2rise(dynFileName,riseFileName,stderr_name)
% dynare2rise(dynFileName,riseFileName,stderr_name,fast)
% dynare2rise(dynFileName,riseFileName,[],fast)
% dynare2rise(dynFileName,riseFileName,stderr_name,[])
% dynare2rise(dynFileName,riseFileName,[],[])

[varargout{1:nargout}]=parser.dynare2rise(varargin{:});

end