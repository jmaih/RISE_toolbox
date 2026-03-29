%  IS_DATE   Check if a value represents a valid date.
% 
%    flag = is_date(x) is a function that checks if the input value (x) represents
%    a valid date. The function utilizes the date2serial function to attempt
%    conversion of the input to a serial date number. If the conversion is successful,
%    the function returns true; otherwise, it returns false.
% 
%    Inputs:
%        - x: The value to be checked for date validity.
% 
%    Output:
%        - flag: A logical value indicating whether the input represents a valid date.
%                It is true if the input is a valid date, and false otherwise.
% 
%    Examples:
%        is_date('2023-06-18')  % Returns true
%        is_date('2023-13-01')  % Returns false
% 
%    Notes:
%        - The function attempts to convert the input value (x) to a serial
%          date number using the date2serial function. If the conversion is
%          successful, it is considered a valid date.
%        - The functions rd, rm, rw, rh, ra, and ry can be used to create valid
%          dates, which can be checked using the is_date function.
%        - For undated dates, simply use integers.
% 
%    See also:
%        date2serial, rd, rm, rw, rh, ra, ry
%