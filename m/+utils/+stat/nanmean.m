% NANMEAN Mean value, ignoring NaNs.
%    M = NANMEAN(X) returns the sample mean of X, treating NaNs as missing
%    values.  For vector input, M is the mean value of the non-NaN elements
%    in X.  For matrix input, M is a row vector containing the mean value of
%    non-NaN elements in each column.  For N-D arrays, NANMEAN operates
%    along the first non-singleton dimension.
% 
%    NANMEAN(X,DIM) takes the mean along dimension DIM of X.
% 
%    See also MEAN, NANMEDIAN, NANSTD, NANVAR, NANMIN, NANMAX, NANSUM.
%
%    Reference page in Doc Center
%       doc stats/nanmean
%
%    Other functions named nanmean
%
%       distributed/nanmean    fints/nanmean
%