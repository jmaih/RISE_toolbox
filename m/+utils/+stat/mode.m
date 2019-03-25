% MODE   Mode, or most frequent value in a sample.
%    M=MODE(X) for vector X computes M as the sample mode, or most frequently
%    occurring value in X.  For a matrix X, M is a row vector containing
%    the mode of each column.  For N-D arrays, MODE(X) is the mode of the
%    elements along the first non-singleton dimension of X.
% 
%    When there are multiple values occurring equally frequently, MODE
%    returns the smallest of those values.  For complex inputs, this is taken
%    to be the first value in a sorted list of values.
% 
%    [M,F]=MODE(X) also returns an array F, of the same size as M.
%    Each element of F is the number of occurrences of the corresponding
%    element of M.
% 
%    [M,F,C]=MODE(X) also returns a cell array C, of the same size
%    as M.  Each element of C is a sorted vector of all the values having
%    the same frequency as the corresponding element of M.
% 
%    [...]=MODE(X,DIM) takes the mode along the dimension DIM of X.
% 
%    This function is most useful with discrete or coarsely rounded data.
%    The mode for a continuous probability distribution is defined as
%    the peak of its density function.  Applying the MODE function to a
%    sample from that distribution is unlikely to provide a good estimate
%    of the peak; it would be better to compute a histogram or density
%    estimate and calculate the peak of that estimate.  Also, the MODE
%    function is not suitable for finding peaks in distributions having
%    multiple modes.
% 
%    Example:
%        X = [3 3 1 4; 0 0 1 1; 0 1 2 4]
%        mode(X,1)
%        mode(X,2)
% 
%       % To find the mode of a continuous variable grouped into bins:
%       y = randn(1000,1);
%       edges = -6:.25:6;
%       bin = discretize(y,edges);
%       m = mode(bin);
%       edges([m, m+1])
%       histogram(y,edges)
% 
%    Class support for input X:
%       float:  double, single
%       integer: uint8, int8, uint16, int16, uint32, int32, uint64, int64
% 
%    See also MEAN, MEDIAN, HISTOGRAM, HISTCOUNTS.
%
%    Reference page in Doc Center
%       doc mode
%
%    Other functions named mode
%
%       categorical/mode      datetime/mode    gpuArray/mode      ts/mode
%       codistributed/mode    duration/mode    timeseries/mode
%