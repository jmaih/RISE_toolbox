%--- help for checksiz ---
%
% CHECKSIZ Test size consistency among input arguments.
% 
%    [ECODE, ERRMSG] = checksiz(SIZES, FUN)
% 
%    Inputs:
%     SIZES - An Mx2 matrix containing the size of each input argument.
% 
%       FUN - The name of the function from which CHECKSIZ is invoked.
% 
%    Outputs:
%     ECODE - Specifies if any size consistencies are violated.
%             1 - Yes
%             0 - No
% 
%    ERRMSG - Displays an error message as an output (instead of to the command
%                  window) if boundary conditions are not met.
% 
%    For example:
%       To check for input argument size consistency in the function FOO which
%       takes the inputs A and B.
% 
%       A = rand(5); B = rand(7);
%       [ecode, errMsg] = checksiz([size(A); size(B)], 'FOO')
% 
%       ecode =
%            1
%  
%       errMsg =
%       Dimensions of inputs are inconsistent.
% 
%    See also CHECKRNG, CHECKTYP.
%