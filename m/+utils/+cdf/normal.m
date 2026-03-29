%--- normal.m not found. Showing help for norm instead. ---
%
% NORM   Vector and matrix norms.
%    N = NORM(V) returns the 2-norm or Euclidean norm of the vector V and is
%    the same as NORM(V,2).
% 
%    N = NORM(V,p) returns the generalized vector p-norm.
%    * If p is a positive real scalar, the p-norm is defined as
%      SUM(ABS(V).^p)^(1/p).
%    * If p = Inf, then N is the largest element of ABS(V).
%    * If p = -Inf, then N is the smallest element of ABS(V).
% 
%    N = NORM(M) returns the 2-norm of the matrix M, which is same as
%    NORM(M,2), and is defined by the largest singular value.
% 
%    N = NORM(M,p) returns the p-norm of the matrix M, where p is 1, 2, or
%    Inf.
%    * If p = 1, then N is the maximum absolute column sum in M.
%    * If p = 2, then N is the maximum singular value of M.
%    * If p = Inf, then N is the maximum absolute row sum of M.
% 
%    N = NORM(X,"fro") returns the Frobenius norm of any numeric array X.
%    Notice that the Frobenius norm for vectors is equivalent to the 2-norm.
%    For N-D arrays, only the Frobenius norm is supported via NORM.
% 
%    By convention, NaN is returned if the first input contains NaNs.
% 
%    See also vecnorm, pagenorm, normest, normalize, cond, hypot.
%
%    Documentation for norm
%       doc norm
%
%