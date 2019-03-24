%--- help for ts.dummy ---
%
%  Creates a dummy observation times series
% 
%  ::
% 
%     db=dummy(start_date,end_date,dummy_date)
% 
%  Args:
% 
%     start_date (char | serial date): start date of the time series
%     end_date (char | serial date): end date of the time series
%     dummy_date (char | serial date): date(s) at which to time series
%       is 1 and not 0
% 
%  Returns:
% 
%     - **db** [ts]: scalar time series of dummy observations
% 
%