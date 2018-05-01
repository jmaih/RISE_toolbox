function db=pageify(pivot_date,varargin)
% pageify -- turn databases into one db with one observation and many pages
%
% ::
%
%
%   db=pageify(pivot_date,varargin)
%
% Args:
%
%    - **pivot_date** [char|serial date] : reference date, typically end of
%    history.
%
%    - **varargin** [struct|ts] : databases in time series format or in struct
%
%    format
%
% Returns:
%    :
%
%    - **db** [ts] : time series with one observation and many pages
%
% Note:
%
%    - This routine is useful for preparing data for conditional forecasting
%
% Example:
%
%    See also:

db=struct2pages(varargin{:});

if db.NumberOfPages>1
    error('Number of pages cannot exceed 1')
end

pivot_date=date2serial(pivot_date);

db=db(pivot_date:db.date_numbers(end));

tmp=double(db);

db=ts(pivot_date,reshape(tmp.',...
    [1,db.NumberOfVariables,db.NumberOfObservations]),...
    db.varnames);

end