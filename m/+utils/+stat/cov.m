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
%    C = cov(...,NANFLAG) specifies how NaN (Not-A-Number) values are 
%    treated. The default is 'includenan':
% 
%    'includenan'   - if the input contains NaN, the output also contains NaN.
%                     Specifically, C(I, J) is NaN if column I or J of X 
%                     contains NaN values.
%    'omitrows'     - omit all rows of X that contain NaN values:
%                       ind = all(~isnan(X), 2);
%                       C = cov(X(ind, :));
%    'partialrows'  - compute each element C(I,J) separately, based only on
%                     the columns I and J of X. Omit rows only if they
%                     contain NaN values in column I or J of X.
%                     The resulting matrix C may not be a positive definite.
%                       ind = all(~isnan(X(:, [I J])));
%                       Clocal = cov(X(ind, [I J]));
%                       C(I, J) = Clocal(1, 2);
% 
%    Class support for inputs X,Y:
%       float: double, single
% 
%    See also CORRCOEF, VAR, STD, MEAN.
%
%    Reference page in Doc Center
%       doc cov
%
%    Other functions named cov
%
%       codistributed/cov    gpuArray/cov    tall/cov    ts/cov
%       fints/cov
%