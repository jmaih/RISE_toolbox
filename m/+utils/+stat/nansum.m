% NANSUM Sum, ignoring NaNs.
%    nansum is not recommended. Use sum instead.
% 
%    Y = NANSUM(X) returns the sum of X, treating NaNs as missing values.
%    For vector input, Y is the sum of the non-NaN elements in X.  For
%    matrix input, Y is a row vector containing the sum of non-NaN elements
%    in each column.  For N-D arrays, NANSUM operates along the first
%    non-singleton dimension.
% 
%    Y = NANSUM(X,'all') sums all of the elements of X.
% 
%    Y = NANSUM(X,DIM) takes the sum along dimension DIM of X.
% 
%    Y = NANSUM(X,VECDIM) sums the elements of X based on the dimensions
%    specified in the vector VECDIM.
% 
%    See also SUM, NANMEAN, NANVAR, NANSTD, NANMIN, NANMAX, NANMEDIAN.
%
%    Documentation for nansum
%       doc nansum
%
%    Other uses of nansum
%
%       distributed/nansum
%