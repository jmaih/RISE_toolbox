% MAX    Largest component.
%    For vectors, MAX(X) is the largest element in X. For matrices,
%    MAX(X) is a row vector containing the maximum element from each
%    column. For N-D arrays, MAX(X) operates along the first
%    non-singleton dimension.
% 
%    [Y,I] = MAX(X) returns the indices of the maximum values in vector I.
%    If the values along the first non-singleton dimension contain more
%    than one maximal element, the index of the first one is returned.
% 
%    MAX(X,Y) returns an array with the largest elements taken from X or Y.
%    X and Y must have compatible sizes. In the simplest cases, they can be
%    the same size or one can be a scalar. Two inputs have compatible sizes
%    if, for every dimension, the dimension sizes of the inputs are either
%    the same or one of them is 1.
% 
%    [Y,I] = MAX(X,[],DIM) operates along the dimension DIM. 
% 
%    When X is complex, the maximum is computed using the magnitude
%    MAX(ABS(X)). In the case of equal magnitude elements, then the phase
%    angle MAX(ANGLE(X)) is used.
% 
%    MAX(...,NANFLAG) specifies how NaN (Not-A-Number) values are treated.
%    NANFLAG can be:
%    'omitnan'    - Ignores all NaN values and returns the maximum of the 
%                   non-NaN elements.  If all elements are NaN, then the
%                   first one is returned.
%    'includenan' - Returns NaN if there is any NaN value.  The index points
%                   to the first NaN element.
%    Default is 'omitnan'.
% 
%    Example: 
%        X = [2 8 4; 7 3 9]
%        max(X,[],1)
%        max(X,[],2)
%        max(X,5)
% 
%    See also MIN, CUMMAX, MEDIAN, MEAN, SORT.
%
%    Reference page in Doc Center
%       doc max
%
%    Other functions named max
%
%       categorical/max      duration/max    sym/max     timeseries/max
%       codistributed/max    fints/max       tall/max    ts/max
%       datetime/max         gpuArray/max
%