function db=nan(start_date,varargin)
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


data=nan(varargin{:});

db=ts(start_date,data,[],[],true);

