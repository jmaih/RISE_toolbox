%  RD DailyDate constructor.
% 
%  Usage:
% 
%    obj = rd(year, month, dayOfTheMonth) creates a
%      rd object for the specified year, month, and day.
% 
%    obj = rd(year, dayOfTheYear) creates a rd
%      object for the specified year and day of the year.
% 
%    obj = rd(dateString) creates a rd object from
%      a string in the format 'mm/dd/yyyy' or 'yyyyDddd' or 'yyyy-mm-dd'.
% 
%    obj = rd(dateNumber) creates a rd object from
%      a serial date number.
% 
%  Inputs:
% 
%    - year: Year of the date.
% 
%    - month: Month of the date.
% 
%    - dayOfTheMonth: Day of the month.
% 
%    - dayOfTheYear: Day of the year.
% 
%    - dateString: A string in the format 'mm/dd/yyyy' or 'yyyyDddd'.
% 
%    - dateNumber: A serial date number representing the date.
% 
%  Example:
% 
%    dateObj = rd(2023, 7, 15);
% 
%    dateObj = rd(2023, 165);
% 
%    dateObj = rd('07/15/2023');
% 
%    dateObj = rd('2023D165');
% 
%    dateObj = rd(737999);
% 
%  See also  rd, rw, rm, rq, rh, ra, ry
%