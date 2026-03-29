%  RH constructor.
% 
%  Usage:
% 
%    obj = rh(year, semi)
% 
%    obj = rh(dateString)
% 
%    obj = rh(datenumber)
% 
%  Inputs:
% 
%    - year: The year of the date.
% 
%    - semi: The semi-annual period (1 or 2).
% 
%    - dateString: A string in the format 'mm/dd/yyyy' or
%      'yyyyHp' where p=1,2 or 'yyyy-mm-dd'.
% 
%    - datenumber: regular date number
% 
%  Example:
% 
%    dateObj = rh(2023, 1);
% 
%    dateObj = rh('2023H1');
% 
%    dateObj = rh('01/01/2023');
% 
%    dateObj = rh(737999);
% 
%  See also rd, rw, rm, rq, ra, ry
%