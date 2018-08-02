function varargout=dynare2rise(varargin)
% Converts dynare model file into rise model files
%
% ::
%
%    dynare2rise(dynFileName);
%    dynare2rise(dynFileName, riseFileName);
%    dynare2rise(dynFileName, riseFileName, stderr_name);
%    dynare2rise(dynFileName, riseFileName, stderr_name, fast);
%
% Args:
%    dynFileName (string): File name of the dynare file to convert
%    riseFileName (string): File name for the output. (defaut: Change file extension to .rs)
%    stderr_name (string): File name to save error messages (default: 'std')
%    fast (bool): (default: true)
%
% Returns:
%    : none
%
% Warning:
%    Dynare has many extensions, and it is hard to define what the core set is. This function converts the basic syntax of a dynare file correctly. However, one should read through the resulting file once to make sure that everything is in order.
%

[varargout{1:nargout}]=parser.dynare2rise(varargin{:});

end