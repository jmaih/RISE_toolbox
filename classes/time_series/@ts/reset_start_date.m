function this=reset_start_date(this,startdate)
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


this=ts(startdate,this.data,this.varnames);
end
