function varargout=mode(this,varargin)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

nout=nargout;
[varargout{1:nout}]=utils.stat.mode(this.data,varargin{:});
end
