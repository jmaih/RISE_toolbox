% NANMAX Maximum value, ignoring NaNs.
%    nanmax is not recommended. Use max instead.
% 
%    M = NANMAX(A) returns the maximum of A with NaNs treated as missing. 
%    For vectors, M is the largest non-NaN element in A.  For matrices, M is
%    a row vector containing the maximum non-NaN element from each column.
%    For N-D arrays, NANMAX operates along the first non-singleton
%    dimension.
% 
%    [M,NDX] = NANMAX(A) returns the indices of the maximum values in A.  If
%    the values along the first non-singleton dimension contain more than
%    one maximal element, the index of the first one is returned.
%   
%    M = NANMAX(A,B) returns an array the same size as A and B with the
%    largest elements taken from A or B.  Either one can be a scalar.
% 
%    M = NANMAX(A,[],'all') returns the largest element of A.
% 
%    [M,NDX] = NANMAX(A,[],DIM) operates along the dimension DIM.
% 
%    M = NANMAX(A,[],VECDIM) returns the largest elements of A based 
%    on the dimensions specified in the vector VECDIM.
% 
%    See also MAX, NANMIN, NANMEAN, NANMEDIAN, NANMIN, NANVAR, NANSTD.
%
%    Documentation for nanmax
%       doc nanmax
%
%    Other uses of nanmax
%
%       distributed/nanmax
%