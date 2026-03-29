%  RW constructor.
% 
%  Usage:
% 
%    obj = rw(year, fullWeeks)
% 
%    obj = rw(dateString)
% 
%    obj = rw(datenumber)
% 
%  Inputs:
% 
%    - year: Year of the date.
% 
%    - fullWeeks: Number of full weeks.
% 
%    - dateString: A string in the format 'mm/dd/yyyy' or 'yyyyWp' or 'yyyy-mm-dd'.
% 
%  Example:
% 
%    dateObj = rw(2023, 10);
% 
%    dateObj = rw('07/01/2023');
% 
%    dateObj = rw('1990W1');
% 
%    dateObj = rw(737999);
% 
% 
%  See also rd, rm, rq, rh, ra, ry
%