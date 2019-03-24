%--- help for ts/ma_filter ---
%
%  Moving average filter
% 
%  ::
% 
%    [trend,detrended]=ma_filter(y,q)
%    [trend,detrended]=ma_filter(y,q,extend)
% 
%  Args:
% 
%     y (ts object): time series object
% 
%     q (integer | {0.5*frequency}): number of periods before or after the
%       current one to be considered in the moving average calculation. The
%       total window length is 2q+1
% 
%     extend (true | {'false'}): if true, replicated observations are
%       added both at the beginning and at the end of the original dataset in
%       order to avoid losing some observations during the filtering process.
% 
%  Returns:
%     :
% 
%     - **trend** [ts] : (non-parametric) trend
%     - **detrended** [ts] : y-trend
% 
%