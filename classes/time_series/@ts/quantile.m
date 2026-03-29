% QUANTILE Quantiles of a sample.
%    Y = QUANTILE(X,P) returns quantiles of the values in X. P is a scalar
%    or a vector of cumulative probability values. When X is a vector, Y is
%    the same size as P, and Y(i) contains the P(i)-th quantile. When X is
%    a matrix, the i-th row of Y contains the P(i)-th quantiles of each
%    column of X. For N-D arrays, QUANTILE operates along the first
%    non-singleton dimension.
% 
%    Y = QUANTILE(X,N) returns quantiles at the N evenly-spaced cumulative
%    probabilities (1:N)/(N+1). N is an integer value greater than 1.
% 
%    Y = QUANTILE(...,"all") calculates quantiles of all the elements in X
%    for either of the first two syntaxes.
% 
%    Y = QUANTILE(...,DIM) calculates quantiles of elements of X along the
%    dimension DIM for either of the first two syntaxes.
% 
%    Y = QUANTILE(...,VECDIM) calculates quantiles of elements of X along
%    the dimensions specified in the vector VECDIM for either of the first
%    two syntaxes. For example, QUANTILE(X,P,[1 2]) operates on the elements
%    contained in the first and second dimensions of X.
% 
%    Y = QUANTILE(...,'Method',METHOD) specifies the method for calculating
%    quantiles. The value of METHOD must be:
%      'exact' - (default) Computes using sorting as explained below. 
%      'approximate' - Computes using an approximation algorithm based on
%      t-digests.
% 
%    Quantiles are specified using cumulative probabilities, from 0 to 1.
%    For an N element vector X, QUANTILE computes quantiles as follows:
%       1) The sorted values in X are taken as the (0.5/N), (1.5/N),
%          ..., ((N-0.5)/N) quantiles.
%       2) Linear interpolation is used to compute quantiles for
%          probabilities between (0.5/N) and ((N-0.5)/N).
%       3) The minimum or maximum values in X are assigned to quantiles
%          for probabilities outside that range.
% 
%    QUANTILE treats NaN values as missing values, and removes them.
% 
%    Examples:
%       y = quantile(x,.50); % the median of x
%       y = quantile(x,[.25 .50 .75]); % the quartiles of x
%       y = quantile(x,3); % another way to get the quartiles of x
%       y = quantile(x,3,2); % the quartiles for each row of x (dim=2)
%       y = quantile(x,[.025 .25 .50 .75 .975]); % a useful summary of x
% 
%    See also IQR, MEDIAN, PRCTILE.
%
%    Documentation for quantile
%       doc quantile
%
%    Other uses of quantile
%
%       ts/quantile
%