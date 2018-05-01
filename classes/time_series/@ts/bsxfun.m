function data=bsxfun(db,fun,b)
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


data=bsxfun(fun,db.data,b);

data=ts(db.start,data);

end