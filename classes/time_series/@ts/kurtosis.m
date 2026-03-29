% KURTOSIS Kurtosis.
%    K = KURTOSIS(X) returns the sample kurtosis of the values in X.  For a
%    vector input, K is the fourth central moment of X, divided by fourth
%    power of its standard deviation.  For a matrix input, K is a row vector
%    containing the sample kurtosis of each column of X.  For N-D arrays,
%    KURTOSIS operates along the first non-singleton dimension.
% 
%    KURTOSIS(X,0) adjusts the kurtosis for bias.  KURTOSIS(X,1) is the same
%    as KURTOSIS(X), and does not adjust for bias.
% 
%    KURTOSIS(X,FLAG,'all') is the kurtosis of all the elements of X.
% 
%    KURTOSIS(X,FLAG,DIM) takes the kurtosis along dimension DIM of X.
% 
%    KURTOSIS(X,FLAG,VECDIM) finds the kurtosis of the elements of X based
%    on the dimensions specified in the vector VECDIM.
% 
%    KURTOSIS treats NaNs as missing values, and removes them.
% 
%    See also MEAN, MOMENT, STD, VAR, SKEWNESS.
%
%    Documentation for kurtosis
%       doc kurtosis
%
%    Other uses of kurtosis
%
%       distributed/kurtosis    tall/kurtosis    ts/kurtosis
%