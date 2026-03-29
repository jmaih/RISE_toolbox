%  DATE2OBS Convert dates to observation numbers based on a specified start date.
% 
%  Inputs:
%  - start_date: The starting date of the time series in RISE format.
%  - new_date: The dates to be converted to observation numbers.
% 
%  Outputs:
%  - obs: The corresponding observation numbers.
%  - flag_: A flag indicating if the new date occurs before the start date.
% 
%  Example:
%    start_date = rq(1990, 1);      % Start date as the first quarter of 1990
%    new_date = rq(1992, 3);        % New date as the third quarter of 1992
%    [obs, flag] = date2obs(start_date, new_date);
% 
%  Note:
%  - RISE assumes that there are no missing observations, and therefore it
%    is enough to provide the first date.
%  - The start_date and new_date arguments should be in RISE format, such as
%    rd(yyyy,mm,dd) for days, rw(yyyy,ww) for weeks, rm(yyyy,mm) for months,
%    rq(yyyy,qq) for quarters, rh(yyyy,hh) for half years, ry(yyyy) for annual
%    or yearly, and integers (double) for undated time series.
%  - The obs output will be NaN for any new date that occurs before the start date.
%  - The flag_ output will be 1 if any new date occurs before the start date,
%    and 0 otherwise.
% 
%  See also: OBS2DATE
%