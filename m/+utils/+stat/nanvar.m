% NANVAR Variance, ignoring NaNs.
%    Y = NANVAR(X) returns the sample variance of the values in X, treating
%    NaNs as missing values.  For a vector input, Y is the variance of the
%    non-NaN elements of X.  For a matrix input, Y is a row vector
%    containing the variance of the non-NaN elements in each column of X.
%    For N-D arrays, NANVAR operates along the first non-singleton dimension
%    of X.
% 
%    NANVAR normalizes Y by N-1 if N>1, where N is the sample size of the 
%    non-NaN elements.  This is an unbiased estimator of the variance of the
%    population from which X is drawn, as long as X consists of independent,
%    identically distributed samples, and data are missing at random.  For
%    N=1, Y is normalized by N. 
% 
%    Y = NANVAR(X,1) normalizes by N and produces the second moment of the
%    sample about its mean.  NANVAR(X,0) is the same as NANVAR(X).
% 
%    Y = NANVAR(X,W) computes the variance using the weight vector W.  The
%    length of W must equal the length of the dimension over which NANVAR
%    operates, and its non-NaN elements must be nonnegative.  Elements of X
%    corresponding to NaN elements of W are ignored.
% 
%    Y = NANVAR(X,W,DIM) takes the variance along dimension DIM of X.
% 
%    See also VAR, NANSTD, NANMEAN, NANMEDIAN, NANMIN, NANMAX, NANSUM.
%
%    Reference page in Doc Center
%       doc stats/nanvar
%
%    Other functions named nanvar
%
%       distributed/nanvar    fints/nanvar
%