%  OBS2DATE Convert observation numbers to dates based on a specified start date.
% 
%  Inputs:
%  - start_date: The starting date of the time series in RISE format.
%  - obs: The observation numbers to be converted to dates.
% 
%  Output:
%  - new_date: The corresponding dates in RISE format.
% 
%  Example:
%    start_date = rq(1990, 1);  % Start date as the first quarter of 1990
%    obs = [1; 2; 3; 4; 5];     % Observation numbers
%    new_date = obs2date(start_date, obs);
% 
%  Note:
%  - RISE assumes that there are no missing observations, and therefore it
%    is enough to provide the first date.
%  - The start_date argument should be in RISE format, such as rd(yyyy,mm,dd)
%    for days, rw(yyyy,ww) for weeks, rm(yyyy,mm) for months, rq(yyyy,qq)
%    for quarters, rh(yyyy,hh) for half years, ry(yyyy) for annual or yearly,
%    and integers (double) for undated time series.
%  - The obs argument should be a positive or negative integer array.
% 
%  See also: DATE2SERIAL
%