% PRCTILE Percentiles of a sample.
%    Y = PRCTILE(X,P) returns percentiles of the values in X. P is a scalar
%    or a vector of percentage values for which to compute percentiles. When
%    X is a vector, Y is the same size as P, and Y(i) contains the P(i)-th
%    percentile. When X is a matrix, the i-th row of Y contains the P(i)-th
%    percentiles of each column of X. For N-D arrays, PRCTILE(X,P) operates
%    along the first non-singleton dimension.
% 
%    Y = PRCTILE(X,P,"all") calculates percentiles of all the elements in X.
% 
%    Y = PRCTILE(X,P,DIM) calculates percentiles along dimension DIM.
% 
%    Y = PRCTILE(X,P,VECDIM) calculates percentiles of elements of X based
%    on the dimensions specified in vector VECDIM. For example,
%    PRCTILE(X,P,[1 2]) operates on the elements contained in the first and
%    second dimensions of X.
% 
%    Y = PRCTILE(...,"Method",METHOD) specifies the method for calculating
%    percentiles. The value of METHOD must be:
%      "exact" - (default) Computes using sorting as explained below.
%      "approximate" - Computes using an approximation algorithm based on
%      t-digests.
% 
%    Percentiles are specified using percentages, from 0 to 100.  For an N
%    element vector X, PRCTILE computes percentiles as follows:
%       1) The sorted values in X are taken as the 100*(0.5/N), 100*(1.5/N),
%          ..., 100*((N-0.5)/N) percentiles.
%       2) Linear interpolation is used to compute percentiles for percent
%          values between 100*(0.5/N) and 100*((N-0.5)/N).
%       3) The minimum or maximum values in X are assigned to percentiles
%          for percent values outside that range.
% 
%    PRCTILE treats NaN values as missing values, and removes them.
% 
%    Examples:
%       y = prctile(x,50); % the median of x
%       y = prctile(x,[2.5 25 50 75 97.5]); % a useful summary of x
% 
%    See also IQR, MEDIAN, QUANTILE.
%
%    Documentation for prctile
%       doc prctile
%
%    Other uses of prctile
%
%       distributed/prctile    tall/prctile    ts/prctile
%