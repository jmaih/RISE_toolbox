function this=reset_start_date(this,startdate)
% Reset the start date
%
% ::
%
%    db = reset_start_date(db,startdate)
%
% Args:
%    db (ts object): original time series object
%    startdate: refer to :func:`ts` to see what startdate formats are available
%
% Returns:
%    :
%    db (ts object): time series object with the start date replaced
%

this=ts(startdate,this.data,this.varnames);
end
