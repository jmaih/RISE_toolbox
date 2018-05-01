function [hms,h,m,s] = sec2hms(secs)
% sec2hms -- converts seconds to hours, minutes and seconds
%
% ::
%
%
%   hms = sec2hms(secs)
%
% Args:
%
%    - **secs** [numeric]:
%
% Returns:
%    :
%
%    - **hms** [char]: string breaking down the input into hours, minutes and
%    seconds.
%
%    - **h** [numeric]:hours.
%
%    - **m** [numeric]: minutes.
%
%    - **s** [numeric]:  seconds.
%
% Note:
%
% Example:
%
%    See also:

% Adapted for RISE from code from Peter J. Acklam 
% http://home.online.no/~pjacklam

secs = round(secs);
s = rem(secs, 60);
m = rem(floor(secs / 60), 60);
h = floor(secs / 3600);
hms = sprintf('%dh%02dm%02ds', h, m, s);
end
