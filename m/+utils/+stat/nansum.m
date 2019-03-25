% NANSUM Sum, ignoring NaNs.
%    Y = NANSUM(X) returns the sum of X, treating NaNs as missing values.
%    For vector input, Y is the sum of the non-NaN elements in X.  For
%    matrix input, Y is a row vector containing the sum of non-NaN elements
%    in each column.  For N-D arrays, NANSUM operates along the first
%    non-singleton dimension.
% 
%    Y = NANSUM(X,DIM) takes the sum along dimension DIM of X.
% 
%    See also SUM, NANMEAN, NANVAR, NANSTD, NANMIN, NANMAX, NANMEDIAN.
%
%    Reference page in Doc Center
%       doc stats/nansum
%
%    Other functions named nansum
%
%       distributed/nansum    fints/nansum
%