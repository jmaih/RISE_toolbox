% MIN    Minimum elements of an array.
%    M = MIN(X) is the smallest element in the vector X. If X is a matrix,
%    M is a row vector containing the minimum element from each column. For 
%    N-D arrays, MIN(X) operates along the first non-singleton dimension.
% 
%    When X is complex, the minimum is computed using the magnitude
%    MIN(ABS(X)). In the case of equal magnitude elements the phase angle 
%    MIN(ANGLE(X)) is used.
% 
%    [M,I] = MIN(X) also returns the indices corresponding to the minimum
%    values. The values in I index into the dimension of X that is being
%    operated on. If X contains more than one element with the minimum
%    value, then the index of the first one is returned.
% 
%    C = MIN(X,Y) returns an array with the smallest elements taken from X 
%    or Y. X and Y must have compatible sizes. In the simplest cases, they 
%    can be the same size or one can be a scalar. Two inputs have compatible 
%    sizes if, for every dimension, the dimension sizes of the inputs are 
%    either the same or one of them is 1.
% 
%    M = MIN(X,[],"all") returns the smallest element of X.
% 
%    [M,I] = MIN(X,[],"all") also returns the linear index into X that
%    corresponds to the minimum value over all elements in X.
% 
%    M = MIN(X,[],DIM) or [M,I] = MIN(X,[],DIM) operates along the 
%    dimension DIM.
% 
%    M = MIN(X,[],VECDIM) operates on the dimensions specified in the vector 
%    VECDIM. For example, MIN(X,[],[1 2]) operates on the elements contained
%    in the first and second dimensions of X.
% 
%    C = MIN(...,NANFLAG) specifies how NaN values are treated:
% 
%    "omitmissing" / "omitnan"       -
%                   (default) Ignores all NaN values and returns the minimum
%                   of the non-NaN elements.
%    "includemissing" / "includenan" -
%                   Returns NaN if there is any NaN value.
% 
%    [M,I] = MIN(X,[],...,"linear") returns the linear index into X that
%    corresponds to the minimum value in X.
% 
%    C = MIN(...,'ComparisonMethod',METHOD) specifies how to compare input
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
%        min(X,[],1)
%        min(X,[],2)
%        min(X,5)
% 
%    See also MAX, BOUNDS, CUMMIN, MEDIAN, MEAN, SORT, MINK, ISLOCALMIN.
%
%    Documentation for min
%       doc min
%
%    Other uses of min
%
%       adolm/min               duration/min            sym/min
%       aplanar/min             gpuArray/min            tabular/min
%       calendarDuration/min    rise_dates.dates/min    tall/min
%       categorical/min         rsymbdiff/min           timeseries/min
%       codistributed/min       splanar/min             ts/min
%       datetime/min
%