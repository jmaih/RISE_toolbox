% NANMIN Minimum value, ignoring NaNs.
%    nanmin is not recommended. Use min instead.
% 
%    M = NANMIN(A) returns the minimum of A with NaNs treated as missing. 
%    For vectors, M is the smallest non-NaN element in A.  For matrices, M
%    is a row vector containing the minimum non-NaN element from each
%    column.  For N-D arrays, NANMIN operates along the first non-singleton
%    dimension.
% 
%    [M,NDX] = NANMIN(A) returns the indices of the minimum values in A.  If
%    the values along the first non-singleton dimension contain more than
%    one minimal element, the index of the first one is returned.
%   
%    M = NANMIN(A,B) returns an array the same size as A and B with the
%    smallest elements taken from A or B.  Either one can be a scalar.
% 
%    M = NANMIN(A,[],'all') returns the smallest element of A.
% 
%    [M,NDX] = NANMIN(A,[],DIM) operates along the dimension DIM.
% 
%    M = NANMIN(A,[],VECDIM) returns the smallest elements of A based 
%    on the dimensions specified in the vector VECDIM.
% 
%    See also MIN, NANMAX, NANMEAN, NANMEDIAN, NANVAR, NANSTD.
%
%    Documentation for nanmin
%       doc nanmin
%
%    Other uses of nanmin
%
%       distributed/nanmin
%