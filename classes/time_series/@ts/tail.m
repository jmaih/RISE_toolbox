%--- help for ts/tail ---
%
%  Returns the last few sample dates of the time series
% 
%  ::
% 
%     db = tail(db);
%     db = tail(db,n);
% 
%  Args:
%     db (ts object): time series object
%     n (integer): number of time steps to show
% 
%  Returns:
%     :
%     - **db** (ts object): time series with the last n-time steps
% 
%  Note:
%     - This is similar to the tail function in stata.
% 
%
%    Other functions named tail
%
%       codistributed/tail    gpuArray/tail    tabular/tail    tall/tail
%