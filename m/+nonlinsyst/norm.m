% NORM   Matrix or vector norm.
%      NORM(X,2) returns the 2-norm of X.
% 
%      NORM(X) is the same as NORM(X,2).
% 
%      NORM(X,1) returns the 1-norm of X.
% 
%      NORM(X,Inf) returns the infinity norm of X.
% 
%      NORM(X,'fro') returns the Frobenius norm of X.
% 
%    In addition, for vectors...
% 
%      NORM(V,P) returns the p-norm of V defined as SUM(ABS(V).^P)^(1/P).
% 
%      NORM(V,Inf) returns the largest element of ABS(V).
% 
%      NORM(V,-Inf) returns the smallest element of ABS(V).
% 
%    By convention, NaN is returned if X or V contains NaNs.
% 
%    See also COND, RCOND, CONDEST, NORMEST, HYPOT.
%
%    Reference page in Doc Center
%       doc norm
%
%    Other functions named norm
%
%       codistributed/norm    gpuArray/norm    sym/norm    tall/norm
%