function db=pageify(pivot_date,varargin)
% pageify -- turn databases into one db with one observation and many pages
%
% Syntax
% -------
% ::
%
%   db=pageify(pivot_date,varargin)
%
% Inputs
% -------
%
% - **pivot_date** [char|serial date] : reference date, typically end of
% history.
%
% - **varargin** [struct|ts] : databases in time series format or in struct
%
% format
%
% Outputs
% --------
%
% - **db** [ts] : time series with one observation and many pages
%
% More About
% ------------
%
% - This routine is useful for preparing data for conditional forecasting
%
% Examples
% ---------
%
% See also: 

db=struct2pages(varargin{:});

if db.NumberOfPages>1
    error('Number of pages cannot exceed 1')
end

if ~ischar(pivot_date)
    pivot_date=serial2date(pivot_date);
end

db=db([pivot_date,':',db.finish]);

tmp=double(db);

db=ts(pivot_date,reshape(tmp.',...
    [1,db.NumberOfVariables,db.NumberOfObservations]),...
    db.varnames);

end