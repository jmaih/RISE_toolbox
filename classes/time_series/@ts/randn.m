function db=randn(start_date,varargin)
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


data=randn(varargin{:});

db=ts(start_date,data);

