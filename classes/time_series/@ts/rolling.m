%--- help for ts/rolling ---
%
%  Applies a function to a rolling window of the time series
% 
%  ::
% 
%     h = rolling(db,func,window);
%     h = rolling(db,func,window,varargin);
% 
%  Args:
% 
%     db (ts object): time series object to get data
% 
%     func: function that will apply to the rolling window of the time
%       series. The function is recursively applied to the data t+(1:window)
%       for t=0,1,2,3,...,T-window, where T is the number of observations
% 
%     window: number of observations in the rolling window
% 
%     varargin: additional arguments to func
% 
%