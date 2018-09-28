function varargout=dynaressfile2rise(varargin)
% converts a basic dynare steady state file into a RISE one
%
% dynaressfile2rise(dynSsfileName,param_names)
% dynaressfile2rise(dynSsfileName,param_names,riseSsfileName)

[varargout{1:nargout}]=parser.dynaressfile2rise(varargin{:});

end