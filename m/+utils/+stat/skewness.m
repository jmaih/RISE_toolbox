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
%    SKEWNESS(X,FLAG,'all') is the skewness of all the elements of X.
% 
%    SKEWNESS(X,FLAG,DIM) takes the skewness along dimension DIM of X.
% 
%    SKEWNESS(X,FLAG,VECDIM) finds the skewness of the elements of X based 
%    on the dimensions specified in the vector VECDIM.
% 
%    SKEWNESS treats NaNs as missing values, and removes them.
% 
%    See also MEAN, MOMENT, STD, VAR, KURTOSIS.
%
%    Documentation for skewness
%       doc skewness
%
%    Other uses of skewness
%
%       distributed/skewness    tall/skewness    ts/skewness
%