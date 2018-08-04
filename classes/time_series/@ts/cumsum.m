function db=cumsum(db,dim)
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

if nargin<2
    dim=1;
end

db=ts(db.date_numbers,cumsum(db.data,dim),db.varnames);
end