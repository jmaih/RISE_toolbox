% MIN    Smallest component.
%    For vectors, MIN(X) is the smallest element in X. For matrices,
%    MIN(X) is a row vector containing the minimum element from each
%    column. For N-D arrays, MIN(X) operates along the first
%    non-singleton dimension.
% 
%    [Y,I] = MIN(X) returns the indices of the minimum values in vector I.
%    If the values along the first non-singleton dimension contain more
%    than one minimal element, the index of the first one is returned.
% 
%    MIN(X,Y) returns an array with the smallest elements taken from X or Y.
%    X and Y must have compatible sizes. In the simplest cases, they can be
%    the same size or one can be a scalar. Two inputs have compatible sizes
%    if, for every dimension, the dimension sizes of the inputs are either
%    the same or one of them is 1.
% 
%    [Y,I] = MIN(X,[],DIM) operates along the dimension DIM.
% 
%    When X is complex, the minimum is computed using the magnitude
%    MIN(ABS(X)). In the case of equal magnitude elements, then the phase
%    angle MIN(ANGLE(X)) is used.
% 
%    MIN(...,NANFLAG) specifies how NaN (Not-A-Number) values are treated.
%    NANFLAG can be:
%    'omitnan'    - Ignores all NaN values and returns the minimum of the 
%                   non-NaN elements.  If all elements are NaN, then the
%                   first one is returned.
%    'includenan' - Returns NaN if there is any NaN value.  The index points
%                   to the first NaN element.
%    Default is 'omitnan'.
% 
%    Example: 
%        X = [2 8 4; 7 3 9]
%        min(X,[],1)
%        min(X,[],2)
%        min(X,5)
% 
%    See also MAX, CUMMIN, MEDIAN, MEAN, SORT.
%
%    Reference page in Doc Center
%       doc min
%
%    Other functions named min
%
%       categorical/min      duration/min    sym/min     timeseries/min
%       codistributed/min    fints/min       tall/min    ts/min
%       datetime/min         gpuArray/min
%