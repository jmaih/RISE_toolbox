function db=expanding(db,func,varargin)
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


window=[];

db=ts_roll_or_expand(db,func,window,varargin{:});

end