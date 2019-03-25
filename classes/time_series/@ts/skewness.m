% SKEWNESS Skewness.
%    S = SKEWNESS(X) returns the sample skewness of the values in X.  For a
%    vector input, S is the third central moment of X, divided by the cube
%    of its standard deviation.  For a matrix input, S is a row vector
%    containing the sample skewness of each column of X.  For N-D arrays,
%    SKEWNESS operates along the first non-singleton dimension.
% 
%    SKEWNESS(X,0) adjusts the skewness for bias.  SKEWNESS(X,1) is the same
%    as SKEWNESS(X), and does not adjust for bias.
% 
%    SKEWNESS(X,FLAG,DIM) takes the skewness along dimension DIM of X.
% 
%    SKEWNESS treats NaNs as missing values, and removes them.
% 
%    See also MEAN, MOMENT, STD, VAR, KURTOSIS.
%
%    Reference page in Doc Center
%       doc skewness
%
%    Other functions named skewness
%
%       distributed/skewness    tall/skewness    ts/skewness
%