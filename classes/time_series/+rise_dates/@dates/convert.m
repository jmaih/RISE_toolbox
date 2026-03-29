%--- help for rise_dates.dates/convert ---
%
%  convert Converts the dates to a different frequency.
% 
%  Usage:
% 
%    d = convert(d, freq)
% 
%  Inputs:
% 
%    - d: The date object.
% 
%    - freq: The target frequency to convert the dates to. Valid values are
%            'D' (daily), 'W' (weekly), 'M' (monthly), 'Q' (quarterly),
%            'H' (half-yearly), and 'Y' (yearly).
% 
%  Outputs:
% 
%    - d: The date object with dates converted to the specified frequency.
% 
%  Example:
% 
%    newDateObj = convert(dateObj, 'M');
% 
%  See also rise_dates.dates.
%
%    Other uses of convert
%
%       matlab.net.http.HeaderField/convert          ts/convert
%       matlab.net.http.io.StringConsumer/convert
%