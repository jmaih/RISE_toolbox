function varargout=dynaressfile2rise(varargin)
% Converts dynare model steady state file into rise steady-state file
%
% ::
%
%    dynaressfile2rise(dynFileName, param_names);
%    dynaressfile2rise(dynFileName, param_names, riseFileName);
%
% Args:
%    dynFileName (string): File name of the dynare file to convert
%    param_names (cell of string): parameter names
%    riseFileName (string): File name for the output. (defaut: Change file extension to .rs)
%
% Returns:
%    : none
%
% Warning:
%    Dynare has many extensions, and it is hard to define what the core set is. This function converts the basic syntax of a dynare file correctly. However, one should read through the resulting file once to make sure that everything is in order.
%

[varargout{1:nargout}]=parser.dynaressfile2rise(varargin{:});

end