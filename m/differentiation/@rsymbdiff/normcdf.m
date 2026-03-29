% NORMCDF Normal cumulative distribution function (cdf).
%    P = NORMCDF(X,MU,SIGMA) returns the cdf of the normal distribution with
%    mean MU and standard deviation SIGMA, evaluated at the values in X.
%    The size of P is the common size of X, MU and SIGMA.  A scalar input
%    functions as a constant matrix of the same size as the other inputs.
% 
%    Default values for MU and SIGMA are 0 and 1, respectively.
% 
%    [P,PLO,PUP] = NORMCDF(X,MU,SIGMA,PCOV,ALPHA) produces confidence bounds
%    for P when the input parameters MU and SIGMA are estimates.  PCOV is a
%    2-by-2 matrix containing the covariance matrix of the estimated parameters.
%    ALPHA has a default value of 0.05, and specifies 100*(1-ALPHA)% confidence
%    bounds.  PLO and PUP are arrays of the same size as P containing the lower
%    and upper confidence bounds.
% 
%    [...] = NORMCDF(...,'upper') computes the upper tail probability of the 
%    normal distribution. This can be used to compute a right-tailed p-value. 
%    To compute a two-tailed p-value, use 2*NORMCDF(-ABS(X),MU,SIGMA).
% 
%    See also ERF, ERFC, NORMFIT, NORMINV, NORMLIKE, NORMPDF, NORMRND, NORMSTAT.
%
%    Documentation for normcdf
%       doc normcdf
%
%    Other uses of normcdf
%
%       adolm/normcdf      rsymbdiff/normcdf    splanar/normcdf
%       aplanar/normcdf
%