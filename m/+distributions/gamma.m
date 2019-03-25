% GAMMA Gamma function.
%    Y = GAMMA(X) evaluates the gamma function for each element of X.
%    X must be real.  The gamma function is defined as:
% 
%       gamma(x) = integral from 0 to inf of t^(x-1) exp(-t) dt.
% 
%    The gamma function interpolates the factorial function.  For
%    integer n, gamma(n+1) = n! (n factorial) = prod(1:n).
% 
%    Class support for input X:
%       float: double, single
% 
%    See also GAMMALN, GAMMAINC, GAMMAINCINV, PSI.
%
%    Reference page in Doc Center
%       doc gamma
%
%    Other functions named gamma
%
%       codistributed/gamma    gpuArray/gamma    sym/gamma
%