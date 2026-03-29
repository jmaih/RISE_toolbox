%  RA constructor for YearlyDate.
% 
%  Usage:
% 
%    obj = ra(year)
% 
%    obj = ra(dateString)
% 
%    dateObj = ra(datenumber);
% 
%  Inputs:
% 
%    - year: Year of the date in double (e.g., 2023)
% 
%    - dateString: A string in the format 'mm/dd/yyyy' or 'yyyy' or 'yyyy-mm-dd'
% 
%    - datenumber: A regular matlab datenumber. Note that the
%      function will treat the input as a date number if it has
%      more than 4 digits. A warning will be issued as well
% 
%  Example:
% 
%    dateObj = ra(2023);
% 
%    dateObj = ra('2023');
% 
%    dateObj = ra('01/01/2023');
% 
%    dateObj = ra(737999);
% 
%  See also rd, rw, rm, rq, rh, ry
%