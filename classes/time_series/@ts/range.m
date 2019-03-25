% RANGE  Sample range.
%    Y = RANGE(X) returns the range of the values in X.  For a vector input,
%    Y is the difference between the maximum and minimum values.  For a
%    matrix input, Y is a vector containing the range for each column.  For
%    N-D arrays, RANGE operates along the first non-singleton dimension.
% 
%    RANGE treats NaNs as missing values, and ignores them.
% 
%    Y = RANGE(X,DIM) operates along the dimension DIM.
% 
%    See also BOUNDS, MIN, MAX, IQR, MAD, STD.
%
%    Reference page in Doc Center
%       doc range
%
%    Other functions named range
%
%       distributed/range    tall/range    ts/range
%