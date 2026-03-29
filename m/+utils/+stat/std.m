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
%    Y = STD(X,0,"all") or Y = STD(X,1,"all") returns the standard deviation
%    of all elements of X. A weight of 0 normalizes by N-1 and a weight of 1 
%    normalizes by N.
% 
%    Y = STD(X,W,DIM) takes the standard deviation along the dimension DIM
%    of X.
% 
%    Y = STD(X,0,VECDIM) or Y = STD(X,1,VECDIM) operates on the dimensions 
%    specified in the vector VECDIM. A weight of 0 normalizes by N-1 and a 
%    weight of 1 normalizes by N. For example, STD(X,0,[1 2]) operates on
%    the elements contained in the first and second dimensions of X.
% 
%    [Y,M] = STD(X,...) also returns the mean M of the values in X that 
%    was used to calculate the standard deviation. If Y is the weighted
%    standard deviation, then M is the weighted mean.
% 
%    The standard deviation is the square root of the variance (VAR).
% 
%    STD(...,NANFLAG) specifies how NaN values are treated:
% 
%    "includemissing" / "includenan" -
%                   (default) The standard deviation of a vector containing
%                   NaN values is also NaN.
%    "omitmissing" / "omitnan"       -
%                   Elements of X or W containing NaN values are ignored.
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
%    Documentation for std
%       doc std
%
%    Other uses of std
%
%       datetime/std                       tabular/std
%       distributed/std                    tall/std
%       duration/std                       timeseries/std
%       prob.LoguniformDistribution/std    ts/std
%       prob.NormalDistribution/std
%