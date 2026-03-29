%  DECIPHER - Interprets error codes returned by RISE
% 
%  Syntax:
%    msgout = decipher(code)
% 
%  Inputs:
%    - code (scalar or vector): Scalar or vector of error codes returned by RISE.
% 
%  Outputs:
%    - msgout (char or cellstr): Explanation of the return code as a char if
%      the input is a scalar or as a cellstr if the input is a vector.
% 
%  Note:
% 
%  Example:
%    msgout = decipher(code);
% 
%  See also:
%