% MAX    Maximum elements of an array.
%    M = MAX(X) is the largest element in the vector X. If X is a matrix, M 
%    is a row vector containing the maximum element from each column. For 
%    N-D arrays, MAX(X) operates along the first non-singleton dimension.
% 
%    When X is complex, the maximum is computed using the magnitude
%    MAX(ABS(X)). In the case of equal magnitude elements the phase angle 
%    MAX(ANGLE(X)) is used.
% 
%    [M,I] = MAX(X) also returns the indices corresponding to the maximum
%    values. The values in I index into the dimension of X that is being
%    operated on. If X contains more than one element with the maximum
%    value, then the index of the first one is returned.
% 
%    C = MAX(X,Y) returns an array with the largest elements taken from X or 
%    Y. X and Y must have compatible sizes. In the simplest cases, they can 
%    be the same size or one can be a scalar. Two inputs have compatible 
%    sizes if, for every dimension, the dimension sizes of the inputs are 
%    either the same or one of them is 1.
% 
%    M = MAX(X,[],"all") returns the largest element of X.
% 
%    [M,I] = MAX(X,[],"all") also returns the linear index into X that 
%    corresponds to the maximum value over all elements in X.
% 
%    M = MAX(X,[],DIM) or [M,I] = MAX(X,[],DIM) operates along the 
%    dimension DIM.
% 
%    M = MAX(X,[],VECDIM) operates on the dimensions specified in the vector 
%    VECDIM. For example, MAX(X,[],[1 2]) operates on the elements contained
%    in the first and second dimensions of X.
% 
%    C = MAX(...,NANFLAG) specifies how NaN values are treated:
% 
%    "omitmissing" / "omitnan"       -
%                   (default) Ignores all NaN values and returns the maximum
%                   of the non-NaN elements.  If all elements are NaN, then 
%                   the first one is returned.
%    "includemissing" / "includenan" -
%                   Returns NaN if there is any NaN value. The index points
%                   to the first NaN element.
% 
%    [M,I] = MAX(X,[],...,'linear') returns the linear index into X that 
%    corresponds to the maximum value in X.
% 
%    C = MAX(...,'ComparisonMethod',METHOD) specifies how to compare input
%    values. The value of METHOD must be:
% 
%        "auto" - (default) Compares real numbers according to "real", and
%                 complex numbers according to "abs".
%        "real" - Compares according to REAL(A). Elements with equal real
%                 parts are then sorted by IMAG(A).
%        "abs"  - Compares according to ABS(A). Elements with equal
%                 magnitudes are then sorted by ANGLE(A).
% 
%    Example: 
%        X = [2 8 4; 7 3 9]
%        max(X,[],1)
%        max(X,[],2)
%        max(X,5)
% 
%    See also MIN, BOUNDS, CUMMAX, MEDIAN, MEAN, SORT, MAXK, ISLOCALMAX.
%
%    Documentation for max
%       doc max
%
%    Other uses of max
%
%       adolm/max               duration/max            sym/max
%       aplanar/max             gpuArray/max            tabular/max
%       calendarDuration/max    rise_dates.dates/max    tall/max
%       categorical/max         rsymbdiff/max           timeseries/max
%       codistributed/max       splanar/max             ts/max
%       datetime/max
%