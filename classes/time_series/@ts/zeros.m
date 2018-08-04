function db=zeros(start_date,varargin)
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


data=zeros(varargin{:});

db=ts(start_date,data);

