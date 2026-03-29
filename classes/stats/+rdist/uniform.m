%--- uniform.m not found. Showing help for isuniform instead. ---
%
%    ISUNIFORM   Check if data is uniformly spaced
%    TF = ISUNIFORM(A) returns a logical scalar that is TRUE if a numeric 
%    vector A is uniformly spaced up to roundoff error. A vector is 
%    uniformly spaced if its elements increase or decrease with a constant,
%    finite step size.
%      
%    [TF,STEP] = ISUNIFORM(A) also returns the step size STEP. If A is
%    uniformly spaced, STEP is a scalar equal to the step size of A. If A is
%    not uniformly spaced, STEP is NaN.
% 
%    Examples:
%        Check if a vector is uniformly spaced:
%        isuniform([2 4 6 8])
% 
%        Check if a vector is uniformly spaced and return the step size:
%        [tf, step] = isuniform([3 6 9]); % step is of type double
%        [tf, step] = isuniform(single([3 6 9])); % step is of type single
%        [tf, step] = isuniform(int32([3 6 9])); % step is of type double
% 
%    See also LINSPACE, COLON, ISREGULAR.
%
%    Documentation for isuniform
%       doc isuniform
%
%