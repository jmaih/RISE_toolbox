%--- igamma1.m not found. Showing help for gamma instead. ---
%
% GAMMA Gamma function.
%    Y = GAMMA(X) evaluates the gamma function for each element of X.
%    X must be real.  The gamma function is defined, for positive x, as:
% 
%       gamma(x) = integral from 0 to inf of t^(x-1) exp(-t) dt.
% 
%    The gamma function interpolates the factorial function. For integer n,
% 
%       gamma(n+1) = factorial(n) = prod(1:n).
% 
%    The domain of the gamma function extends to negative real numbers by
%    analytic continuation, with simple poles at the negative integers.
% 
%    Class support for input X:
%       float: double, single
% 
%    See also FACTORIAL, GAMMALN, GAMMAINC, GAMMAINCINV, PSI.
%
%    Documentation for gamma
%       doc gamma
%
%