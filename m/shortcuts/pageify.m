function db=pageify(varargin)
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

db=utils.time_series.pageify(varargin{:});

end