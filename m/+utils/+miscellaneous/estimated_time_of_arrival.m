function [hrs,mins,secs]=estimated_time_of_arrival(start_clock,percentage_complete)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

if percentage_complete<=0||percentage_complete>1
    error('second argument must be in (0,1]')
end
secs = etime(clock,start_clock)*(1-percentage_complete)/percentage_complete;
mins = floor(secs/60);
secs = ceil(secs - 60*mins);
hrs  = floor(mins/60);
mins = mins - hrs*60;
end