%  Highlights or shades specific periods between start and finish years.
% 
%  Syntax:
% <<
%    shade(start_finish)
% 
%    shade(start_finish, colorstr)
% 
%    shade(start_finish, colorstr, fig)
% >>
%  Inputs:
% 
%    - `start_finish`:
%      - Time-series (`ts`) of data with zeros and ones, where the first "1"
%        in a sequence represents the beginning of an area to shade and the
%        last "1" in the same sequence, the end of the area to shade.
%      -  Row represents the start and finish dates of a period. The dates
%        should be in a format recognizable by RISE. Valid formats include
%        integers or members of the `rise_dates.dates` class.
% 
%    - `colorstr`: String, character, or 1x3 RGB vector specifying the fill
%      color. If a string or character is provided, it should be one of the
%      following: 'r', 'g', 'b', 'c', 'm', 'y', 'w', 'k'.  If an RGB vector
%      is provided, it should be a 1x3 vector specifying the color values
%      divided by 255. Default is [211, 211, 211]/255 (light gray).
% 
%    - `fig`: Figure handle or array of figure handles. Default is the
%      current figure. 
% 
%  Outputs:
% 
%    - None
% 
%  Details:
%    The function generates shaded areas in a plot to highlight specific
%    periods between start and finish years. It can be used to visually
%    represent recessions or other specific time intervals.
% 
%  Example usage:
% <<
%    start_finish = [2000, 2002; 2005, 2008; 2010, 2015];
% 
%    shade(start_finish, 'y', gcf);
% >>
%