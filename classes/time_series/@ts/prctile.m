% PRCTILE Percentiles of a sample.
%    Y = PRCTILE(X,P) returns percentiles of the values in X.  P is a scalar
%    or a vector of percent values.  When X is a vector, Y is the same size
%    as P, and Y(i) contains the P(i)-th percentile.  When X is a matrix,
%    the i-th row of Y contains the P(i)-th percentiles of each column of X.
%    For N-D arrays, PRCTILE operates along the first non-singleton
%    dimension.
% 
%    Y = PRCTILE(X,P,DIM) calculates percentiles along dimension DIM.  The
%    DIM'th dimension of Y has length LENGTH(P).
% 
%    Percentiles are specified using percentages, from 0 to 100.  For an N
%    element vector X, PRCTILE computes percentiles as follows:
%       1) The sorted values in X are taken as the 100*(0.5/N), 100*(1.5/N),
%          ..., 100*((N-0.5)/N) percentiles.
%       2) Linear interpolation is used to compute percentiles for percent
%          values between 100*(0.5/N) and 100*((N-0.5)/N)
%       3) The minimum or maximum values in X are assigned to percentiles
%          for percent values outside that range.
% 
%    PRCTILE treats NaNs as missing values, and removes them.
% 
%    Examples:
%       y = prctile(x,50); % the median of x
%       y = prctile(x,[2.5 25 50 75 97.5]); % a useful summary of x
% 
%    See also IQR, MEDIAN, NANMEDIAN, QUANTILE.
%
%    Reference page in Doc Center
%       doc prctile
%
%    Other functions named prctile
%
%       distributed/prctile    ts/prctile
%