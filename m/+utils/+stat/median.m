% MEDIAN Median value.
%    For vectors, MEDIAN(x) is the median value of the elements in x.
%    For matrices, MEDIAN(X) is a row vector containing the median value
%    of each column.  For N-D arrays, MEDIAN(X) is the median value of the
%    elements along the first non-singleton dimension of X.
% 
%    MEDIAN(X,"all") is the median of all elements of X.
% 
%    MEDIAN(X,DIM) takes the median along the dimension DIM of X.
% 
%    MEDIAN(X,VECDIM) operates on the dimensions specified in the vector 
%    VECDIM. For example, MEDIAN(X,[1 2]) operates on the elements contained
%    in the first and second dimensions of X.
% 
%    MEDIAN(...,NANFLAG) specifies how NaN values are treated:
% 
%    "includemissing" / "includenan" -
%                   (default) The median of a vector containing NaN values
%                   is also NaN.
%    "omitmissing" / "omitnan"       -
%                   The median of a vector containing NaN values is the
%                   median of all its non-NaN elements. If all elements
%                   are NaN, the result is NaN.
% 
%    Example:
%        X = [1 2 4 4; 3 4 6 6; 5 6 8 8; 5 6 8 8]
%        median(X,1)
%        median(X,2)
% 
%    Class support for input X:
%       float: double, single
%       integer: uint8, int8, uint16, int16, uint32, int32, uint64, int64
% 
%    See also MEAN, STD, MIN, MAX, VAR, COV, MODE.
%
%    Documentation for median
%       doc median
%
%    Other uses of median
%
%       categorical/median
%       codistributed/median
%       datetime/median
%       duration/median
%       prob.LoguniformDistribution/median
%       prob.NormalDistribution/median
%       tabular/median
%       tall/median
%       timeseries/median
%       ts/median
%