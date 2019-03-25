% STD Standard deviation.
%    For vectors, Y = STD(X) returns the standard deviation.  For matrices,
%    Y is a row vector containing the standard deviation of each column.  For
%    N-D arrays, STD operates along the first non-singleton dimension of X.
% 
%    STD normalizes Y by N-1 if N>1, where N is the sample size.  This is
%    the sqrt of an unbiased estimator of the variance of the population
%    from which X is drawn, as long as X consists of independent,
%    identically distributed samples. For N=1, Y is normalized by N.
% 
%    Y = STD(X,1) normalizes by N and produces the square root of the second
%    moment of the sample about its mean.  STD(X,0) is the same as STD(X).
% 
%    Y = STD(X,W) computes the standard deviation using the weight vector W.
%    W typically contains either counts or inverse variances.  The length of
%    W must equal the length of the dimension over which STD operates, and
%    its elements must be nonnegative.  If X(I) is assumed to have standard
%    deviation proportional to 1/SQRT(W(I)), then Y * SQRT(MEAN(W)/W(I)) is
%    an estimate of the standard deviation of X(I).  In other words, Y *
%    SQRT(MEAN(W)) is an estimate of standard deviation for an observation
%    given weight 1.
% 
%    Y = STD(X,W,DIM) takes the standard deviation along the dimension DIM
%    of X.  Pass in 0 for W to use the default normalization by N-1, or 1 to
%    use N.
% 
%    The standard deviation is the square root of the variance (VAR).
% 
%    STD(...,NANFLAG) specifies how NaN (Not-A-Number) values are treated.
%    The default is 'includenan':
% 
%    'includenan' - the standard deviation of a vector containing NaN 
%                   values is also NaN.
%    'omitnan'    - elements of X or W containing NaN values are ignored.
%                   If all elements are NaN, the result is NaN.
% 
%    Example: 
%        X = [4 -2 1; 9 5 7]
%        std(X,0,1)
%        std(X,0,2)
% 
%    Class support for inputs X, W:
%       float: double, single
% 
%    See also COV, MEAN, VAR, MEDIAN, CORRCOEF.
%
%    Reference page in Doc Center
%       doc std
%
%    Other functions named std
%
%       datetime/std       duration/std    tall/std          ts/std
%       distributed/std    fints/std       timeseries/std
%