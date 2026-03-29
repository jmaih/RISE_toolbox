% NORMINV Inverse of the normal cumulative distribution function (cdf).
%    X = NORMINV(P,MU,SIGMA) returns the inverse cdf for the normal
%    distribution with mean MU and standard deviation SIGMA, evaluated at
%    the values in P.  The size of X is the common size of the input
%    arguments.  A scalar input functions as a constant matrix of the same
%    size as the other inputs.
% 
%    Default values for MU and SIGMA are 0 and 1, respectively.
% 
%    [X,XLO,XUP] = NORMINV(P,MU,SIGMA,PCOV,ALPHA) produces confidence bounds
%    for X when the input parameters MU and SIGMA are estimates.  PCOV is a
%    2-by-2 matrix containing the covariance matrix of the estimated parameters.
%    ALPHA has a default value of 0.05, and specifies 100*(1-ALPHA)% confidence
%    bounds.  XLO and XUP are arrays of the same size as X containing the lower
%    and upper confidence bounds.
% 
%    See also ERFINV, ERFCINV, NORMCDF, NORMFIT, NORMLIKE, NORMPDF,
%             NORMRND, NORMSTAT.
%
%    Documentation for norminv
%       doc norminv
%
%    Other uses of norminv
%
%       distributed/norminv    rsymbdiff/norminv    splanar/norminv
%