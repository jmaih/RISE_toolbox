% COV Covariance matrix.
%    COV(X), if X is a vector, returns the variance.  For matrices, where 
%    each row is an observation, and each column a variable, COV(X) is the 
%    covariance matrix.  DIAG(COV(X)) is a vector of variances for each 
%    column, and SQRT(DIAG(COV(X))) is a vector of standard deviations. 
%    COV(X,Y), where X and Y are matrices with the same number of elements,
%    is equivalent to COV([X(:) Y(:)]). 
%    
%    COV(X) or COV(X,Y) normalizes by (N-1) if N>1, where N is the number of
%    observations.  This makes COV(X) the best unbiased estimate of the
%    covariance matrix if the observations are from a normal distribution.
%    For N=1, COV normalizes by N.
% 
%    COV(X,1) or COV(X,Y,1) normalizes by N and produces the second
%    moment matrix of the observations about their mean.  COV(X,Y,0) is
%    the same as COV(X,Y) and COV(X,0) is the same as COV(X).
% 
%    The mean is removed from each column before calculating the result.
% 
%    C = COV(...,NANFLAG) specifies how NaN values are treated:
% 
%    "includemissing" / "includenan" -
%                     (default) If the input contains NaN, the output also
%                     contains NaN. Specifically, C(I, J) is NaN if column I
%                     or J of X contains NaN values.
%    "omitrows"     - Omits all rows of X that contain NaN values and
%                     returns the variance of all other rows.
%    "partialrows"  - Computes each element C(I,J) separately, based only on
%                     the columns I and J of X. Omits rows only if they
%                     contain NaN values in column I or J of X. The resulting
%                     matrix C may not be positive definite.
% 
%    Class support for inputs X,Y:
%       float: double, single
% 
%    See also CORRCOEF, VAR, STD, MEAN.
%
%    Documentation for cov
%       doc cov
%
%    Other uses of cov
%
%       codistributed/cov    symfun/cov    tall/cov    ts/cov
%       gpuArray/cov
%