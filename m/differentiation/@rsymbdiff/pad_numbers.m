%--- help for rsymbdiff.pad_numbers ---
%
%  PAD_NUMBERS - Convert a numeric array to a formatted string
% 
%    [c, n] = pad_numbers(d) converts a numeric array 'd' into a formatted
%    string 'c' with options for integer and floating-point precision.
%    Additionally, it returns the number of elements in the input array.
% 
%  Input:
%    - d (numeric array): Input numeric array to be converted.
% 
%  Output:
%    - c (char): Formatted string containing the numeric elements.
%    - n (integer): Number of elements in the input array 'd'.
% 
%  Example:
%    d = [1, 2.5, 3.14159];
%    [c, n] = pad_numbers(d);
% 
%  See also: sprintf
%