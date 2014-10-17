function varargout=mode(this,varargin)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

nout=nargout;
[varargout{1:nout}]=utils.stat.mode(this.data,varargin{:});
end
