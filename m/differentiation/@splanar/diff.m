% DIFF Difference and approximate derivative.
%    DIFF(X), for a vector X, is [X(2)-X(1)  X(3)-X(2) ... X(n)-X(n-1)].
%    DIFF(X), for a matrix X, is the matrix of row differences,
%       [X(2:n,:) - X(1:n-1,:)].
%    DIFF(X), for an N-D array X, is the difference along the first
%       non-singleton dimension of X.
%    DIFF(X,N) is the N-th order difference along the first non-singleton 
%       dimension (denote it by DIM). If N >= size(X,DIM), DIFF takes 
%       successive differences along the next non-singleton dimension.
%    DIFF(X,N,DIM) is the Nth difference function along dimension DIM. 
%       If N >= size(X,DIM), DIFF returns an empty array.
% 
%    Examples:
%       h = .001; x = 0:h:pi;
%       diff(sin(x.^2))/h is an approximation to 2*cos(x.^2).*x
%       diff((1:10).^2) is 3:2:19
% 
%       If X = [3 7 5
%               0 9 2]
%       then diff(X,1,1) is [-3 2 -3], diff(X,1,2) is [4 -2
%                                                      9 -7],
%       diff(X,2,2) is the 2nd order difference along the dimension 2, and
%       diff(X,3,2) is the empty matrix.
% 
%    See also GRADIENT, SUM, PROD.
%
%    Reference page in Doc Center
%       doc diff
%
%    Other functions named diff
%
%       codistributed/diff    duration/diff    gpuArray/diff    sym/diff
%       datetime/diff         fints/diff       splanar/diff     tall/diff
%