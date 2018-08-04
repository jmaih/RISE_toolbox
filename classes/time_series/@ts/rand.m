function db=rand(start_date,varargin)
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


data=rand(varargin{:});

db=ts(start_date,data);

