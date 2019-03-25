% MEAN   Average or mean value.
%    S = MEAN(X) is the mean value of the elements in X if X is a vector. 
%    For matrices, S is a row vector containing the mean value of each 
%    column. 
%    For N-D arrays, S is the mean value of the elements along the first 
%    array dimension whose size does not equal 1.
% 
%    MEAN(X,DIM) takes the mean along the dimension DIM of X.
% 
%    S = MEAN(...,TYPE) specifies the type in which the mean is performed, 
%    and the type of S. Available options are:
% 
%    'double'    -  S has class double for any input X
%    'native'    -  S has the same class as X
%    'default'   -  If X is floating point, that is double or single,
%                   S has the same class as X. If X is not floating point, 
%                   S has class double.
% 
%    S = MEAN(...,NANFLAG) specifies how NaN (Not-A-Number) values are 
%    treated. The default is 'includenan':
% 
%    'includenan' - the mean of a vector containing NaN values is also NaN.
%    'omitnan'    - the mean of a vector containing NaN values is the mean 
%                   of all its non-NaN elements. If all elements are NaN,
%                   the result is NaN.
% 
%    Example:
%        X = [1 2 3; 3 3 6; 4 6 8; 4 7 7]
%        mean(X,1)
%        mean(X,2)
% 
%    Class support for input X:
%       float: double, single
%       integer: uint8, int8, uint16, int16, uint32,
%                int32, uint64, int64
% 
%    See also MEDIAN, STD, MIN, MAX, VAR, COV, MODE.
%
%    Reference page in Doc Center
%       doc mean
%
%    Other functions named mean
%
%       codistributed/mean    fints/mean       timeseries/mean
%       datetime/mean         gpuArray/mean    ts/mean
%       duration/mean         tall/mean
%