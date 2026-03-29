%  RQ constructor.
% 
%  Usage:
% 
%    - obj = rq(year, quarter)
% 
%    - obj = rq(dateString)
% 
%    - obj = rq(datenumber)
% 
%  Inputs:
% 
%    - year: Year of the date.
% 
%    - quarter: Quarter of the date (1, 2, 3, or 4).
% 
%    - dateString: A string in the format 'mm/dd/yyyy' or
%      'yyyyQp' where p=1,2,3,4 or 'yyyy-mm-dd'
% 
%  Example:
% 
%    dateObj = rq(2023, 2);
% 
%    dateObj = rq('2023Q2');
% 
%    dateObj = rq('04/01/2023');
% 
%    dateObj = rq(737999);
% 
% 
%  See also rd, rw, rm, rh, ra, ry
%