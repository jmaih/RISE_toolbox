%  pageify -- turn databases into one db with one observation and many pages
% 
%  ::
% 
% 
%    db=pageify(pivot_date,varargin)
% 
%  Args:
% 
%     - **pivot_date** [char|serial date] : reference date, typically end of
%     history.
% 
%     - **varargin** [struct|ts] : databases in time series format or in struct
% 
%     format
% 
%  Returns:
%     :
% 
%     - **db** [ts] : time series with one observation and many pages
% 
%  Note:
% 
%     - This routine is useful for preparing data for conditional forecasting
% 
%  Example:
% 
%     See also:
%