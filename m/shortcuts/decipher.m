function varargout=decipher(varargin)
% decipher - interprets error codes
%
% Syntax
% -------
% ::
%
%   msgout=decipher(code)
%
% Inputs
% -------
%
% - **code** [scalar|vector]: scalar or vector of error codes returned by
% RISE.
%
% Outputs
% --------
%
% - **msgout** [char|cellstr]: explanation of the return code as a char if
% the input is a scalar or as a cellstr if input is a vector.
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

[varargout{1:nargout}]=utils.error.decipher(varargin{:});