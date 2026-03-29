%--- rdist.m not found. Showing help for pdist instead. ---
%
% PDIST Pairwise distance between observations.
%    D = PDIST(X) returns a vector D containing the Euclidean distances
%    between each pair of observations in the M-by-N data matrix X. Rows of
%    X correspond to observations, columns correspond to variables. D is a
%    1-by-(M*(M-1)/2) row vector, corresponding to the M*(M-1)/2 pairs of
%    observations in X.
% 
%    D = PDIST(X, DISTANCE) computes D using DISTANCE.  Choices are:
% 
%        'euclidean'        - Euclidean distance (default)
%        'squaredeuclidean' - Squared Euclidean distance 
%        'seuclidean'       - Standardized Euclidean distance. Each
%                             coordinate difference between rows in X is
%                             scaled by dividing by the corresponding
%                             element of the standard deviation S=NANSTD(X).
%                             To specify another value for S, use
%                             D=PDIST(X,'seuclidean',S).
%        'fasteuclidean'    - Euclidean distance computed by using an
%                             alternative algorithm that saves time. This
%                             faster algorithm can, in some cases, reduce
%                             accuracy.
%        'fastsquaredeuclidean'
%                           - Squared Euclidean distance computed by using 
%                             an alternative algorithm that saves time. This
%                             faster algorithm can, in some cases, reduce
%                             accuracy.
%        'fastseuclidean'   - Standardized Euclidean distance computed by 
%                             using an alternative algorithm that saves 
%                             time. This faster algorithm can, in some 
%                             cases, reduce accuracy.
%        'cityblock'        - City Block distance
%        'minkowski'        - Minkowski distance. The default exponent is 2.
%                             To specify a different exponent, use
%                             D = PDIST(X,'minkowski',P), where the exponent
%                             P is a scalar positive value.
%        'chebychev'        - Chebychev distance (maximum coordinate
%                             difference)
%        'mahalanobis'      - Mahalanobis distance, using the sample
%                             covariance of X as computed by NANCOV. To
%                             compute the distance with a different
%                             covariance, use
%                             D =  PDIST(X,'mahalanobis',C), where the
%                             matrix C is symmetric and positive definite.
%        'cosine'           - One minus the cosine of the included angle
%                             between observations (treated as vectors)
%        'correlation'      - One minus the sample linear correlation
%                             between observations (treated as sequences of
%                             values).
%        'spearman'         - One minus the sample Spearman's rank 
%                             correlation between observations (treated as
%                             sequences of values).
%        'hamming'          - Hamming distance, percentage of coordinates
%                             that differ
%        'jaccard'          - One minus the Jaccard coefficient, the
%                             percentage of nonzero coordinates that differ
%        function           - A distance function specified using @, for
%                             example @DISTFUN.
% 
%    A distance function must be of the form
% 
%          function D2 = DISTFUN(XI, XJ),
% 
%    taking as arguments a 1-by-N vector XI containing a single row of X, an
%    M2-by-N matrix XJ containing multiple rows of X, and returning an
%    M2-by-1 vector of distances D2, whose Jth element is the distance
%    between the observations XI and XJ(J,:).
% 
%    The output D is arranged in the order of ((2,1),(3,1),..., (M,1),
%    (3,2),...(M,2),.....(M,M-1)), i.e. the lower left triangle of the full
%    M-by-M distance matrix in column order.  To get the distance between
%    the Ith and Jth observations (I < J), either use the formula
%    D((I-1)*(M-I/2)+J-I), or use the helper function Z = SQUAREFORM(D),
%    which returns an M-by-M square symmetric matrix, with the (I,J) entry
%    equal to distance between observation I and observation J.
% 
%    D = PDIST(X,DISTANCE,'CacheSize',CACHESIZE) uses an intermediate matrix 
%    stored in cache to compute D, when 'Distance' is one of 
%    {'fasteuclidean','fastsquaredeuclidean','fastseuclidean'}. 'CacheSize'
%    can be a positive scalar or 'maximal'. The default is 1e3.
%    If numeric, 'CacheSize' specifies the cache size in megabytes (MB) to
%    allocate for an intermediate matrix.
%    If 'maximal', pdist attempts to allocate enough memory for an entire
%    intermediate matrix whose size is M-by-M (M is the number of rows of
%    the input data X).
%    'CacheSize' does not have to be large enough for an entire intermediate
%    matrix, but it must be at least large enough to hold an M-by-1 vector. 
%    Otherwise, the regular algorithm of computing Euclidean distance will
%    be used instead. If the specified cache size exceeds the available
%    memory, MATLAB issues an out-of-memory error.
% 
%    Example:
%       % Compute the ordinary Euclidean distance
%       X = randn(100, 5);                 % some random points
%       D = pdist(X, 'euclidean');         % euclidean distance
% 
%       % Compute the Euclidean distance with each coordinate difference
%       % scaled by the standard deviation
%       Dstd = pdist(X,'seuclidean');
% 
%       % Use a function handle to compute a distance that weights each
%       % coordinate contribution differently
%       Wgts = [.1 .3 .3 .2 .1];           % coordinate weights
%       weuc = @(XI,XJ,W)(sqrt((XI-XJ).^2 * W'));
%       Dwgt = pdist(X, @(Xi,Xj) weuc(Xi,Xj,Wgts));
% 
%    See also SQUAREFORM, LINKAGE, SILHOUETTE, PDIST2.
%
%    Documentation for pdist
%       doc pdist
%
%