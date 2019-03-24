%--- help for ts.step_dummy ---
%
%  Implements a step dummy in the time series
% 
%  Args:
% 
%     start_date: start date of the time series
% 
%     end_date: end date of the time series
% 
%     dummy_start_date: start date of the step dummy
% 
%  Returns:
%     :
% 
%     - **db** (ts object): a time series with zeros from the first observation to the date
%       before the start of the dummy and ones from then on.
% 
%