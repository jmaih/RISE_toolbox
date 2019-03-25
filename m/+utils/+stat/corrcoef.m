% CORRCOEF Correlation coefficients.
%    R=CORRCOEF(X) calculates a matrix R of correlation coefficients for
%    an array X, in which each row is an observation and each column is a
%    variable.
% 
%    R=CORRCOEF(X,Y), where X and Y are column vectors, is the same as
%    R=CORRCOEF([X Y]).  CORRCOEF converts X and Y to column vectors if they
%    are not, i.e., R=CORRCOEF(X,Y) is equivalent to R=CORRCOEF([X(:) Y(:)])
%    in that case.
% 
%    If C is the covariance matrix, C = COV(X), then CORRCOEF(X) is
%    the matrix whose (i,j)'th element is
% 
%           C(i,j)/SQRT(C(i,i)*C(j,j)).
% 
%    [R,P]=CORRCOEF(...) also returns P, a matrix of p-values for testing
%    the hypothesis of no correlation.  Each p-value is the probability
%    of getting a correlation as large as the observed value by random
%    chance, when the true correlation is zero.  If P(i,j) is small, say
%    less than 0.05, then the correlation R(i,j) is significant.
% 
%    [R,P,RLO,RUP]=CORRCOEF(...) also returns matrices RLO and RUP, of
%    the same size as R, containing lower and upper bounds for a 95%
%    confidence interval for each coefficient.
% 
%    [...]=CORRCOEF(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies additional
%    parameters and their values.  Valid parameters are the following:
% 
%        Parameter  Value
%        'Alpha'    A number between 0 and 1 to specify a confidence
%                   level of 100*(1-ALPHA)%.  Default is 0.05 for 95%
%                   confidence intervals.
%        'Rows'     Either 'all' (default) to use all rows, 'complete' to
%                   use rows with no NaN values, or 'pairwise' to compute
%                   R(i,j) using rows with no NaN values in column i or j.
%                   The 'pairwise' option potentially uses different sets
%                   of rows to compute different elements of R, and can
%                   produce a matrix that is indefinite.
% 
%    The p-value is computed by transforming the correlation to create a t
%    statistic having N-2 degrees of freedom, where N is the number of rows
%    of X.  The confidence bounds are based on an asymptotic normal
%    distribution of 0.5*log((1+R)/(1-R)), with an approximate variance equal
%    to 1/(N-3).  These bounds are accurate for large samples when X has a
%    multivariate normal distribution.
% 
%    Example:  Generate random data having correlation between column 4
%              and the other columns.
%        x = randn(30,4);       % uncorrelated data
%        x(:,4) = sum(x,2);     % introduce correlation
%        [r,p] = corrcoef(x)    % compute sample correlation and p-values
%        [i,j] = find(p<0.05);  % find significant correlations
%        [i,j]                  % display their (row,col) indices
% 
%    Class support for inputs X,Y:
%       float: double, single
% 
%    See also COV, VAR, STD, PLOTMATRIX.
%
%    Reference page in Doc Center
%       doc corrcoef
%
%    Other functions named corrcoef
%
%       codistributed/corrcoef    gpuArray/corrcoef    ts/corrcoef
%       fints/corrcoef
%