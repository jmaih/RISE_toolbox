% SUM Sum of elements.
%    S = SUM(X) is the sum of the elements of the vector X. If X is a matrix,
%    S is a row vector with the sum over each column. For N-D arrays, 
%    SUM(X) operates along the first non-singleton dimension.
% 
%    S = SUM(X,DIM) sums along the dimension DIM. 
% 
%    S = SUM(...,TYPE) specifies the type in which the 
%    sum is performed, and the type of S. Available options are:
% 
%    'double'    -  S has class double for any input X
%    'native'    -  S has the same class as X
%    'default'   -  If X is floating point, that is double or single,
%                   S has the same class as X. If X is not floating point, 
%                   S has class double.
% 
%    S = SUM(...,NANFLAG) specifies how NaN (Not-A-Number) values are 
%    treated. The default is 'includenan':
% 
%    'includenan' - the sum of a vector containing NaN values is also NaN.
%    'omitnan'    - the sum of a vector containing NaN values
%                   is the sum of all its non-NaN elements. If all 
%                   elements are NaN, the result is 0.
% 
%    Examples:
%        X = [0 1 2; 3 4 5]
%        sum(X, 1)
%        sum(X, 2)
% 
%        X = int8(1:20)
%        sum(X)             % returns double(210), accumulates in double
%        sum(X,'native')    % returns int8(127), because it accumulates in
%                           % int8 but overflows and saturates.
% 
%    See also PROD, CUMSUM, DIFF, ACCUMARRAY, ISFLOAT.
%
%    Reference page in Doc Center
%       doc sum
%
%    Other functions named sum
%
%       codistributed/sum    gpuArray/sum    tall/sum          ts/sum
%       duration/sum         sym/sum         timeseries/sum
%