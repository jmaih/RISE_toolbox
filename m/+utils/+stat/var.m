% VAR Variance.
%    For vectors, Y = VAR(X) returns the variance of the values in X.  For
%    matrices, Y is a row vector containing the variance of each column of
%    X.  For N-D arrays, VAR operates along the first non-singleton
%    dimension of X.
% 
%    VAR normalizes Y by N-1 if N>1, where N is the sample size.  This is
%    an unbiased estimator of the variance of the population from which X is
%    drawn, as long as X consists of independent, identically distributed
%    samples. For N=1, Y is normalized by N. 
% 
%    Y = VAR(X,1) normalizes by N and produces the second moment of the
%    sample about its mean.  VAR(X,0) is the same as VAR(X).
% 
%    Y = VAR(X,W) computes the variance using the weight vector W.  W 
%    typically contains either counts or inverse variances.  The length of W 
%    must equal the length of the dimension over which VAR operates, and its
%    elements must be nonnegative.  If X(I) is assumed to have variance 
%    proportional to 1/W(I), then Y * MEAN(W)/W(I) is an estimate of the 
%    variance of X(I).  In other words, Y * MEAN(W) is an estimate of 
%    variance for an observation given weight 1.
% 
%    Y = VAR(X,W,DIM) takes the variance along the dimension DIM of X.  Pass
%    in 0 for W to use the default normalization by N-1, or 1 to use N.
% 
%    The variance is the square of the standard deviation (STD).
% 
%    VAR(...,NANFLAG) specifies how NaN (Not-A-Number) values are treated.
%    The default is 'includenan':
% 
%    'includenan' - the variance of a vector containing NaN values 
%                   is also NaN.
%    'omitnan'    - elements of X or W containing NaN values are ignored.
%                   If all elements are NaN, the result is NaN.
% 
%    Example:
%        X = [4 -2 1; 9 5 7]
%        var(X,0,1)
%        var(X,0,2)
% 
%    Class support for inputs X, W:
%       float: double, single
% 
%    See also MEAN, STD, COV, CORRCOEF.
%
%    Reference page in Doc Center
%       doc var
%
%    Other functions named var
%
%       codistributed/var    gpuArray/var    timeseries/var    ts/var
%       fints/var            tall/var
%