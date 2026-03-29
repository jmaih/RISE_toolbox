% BETAINV Inverse of the beta cumulative distribution function (cdf).
%    X = BETAINV(P,A,B) returns the inverse of the beta cdf with 
%    parameters A and B at the values in P.
% 
%    The size of X is the common size of the input arguments. A scalar input  
%    functions as a constant matrix of the same size as the other inputs.    
% 
%    BETAINV uses Newton's method to converge to the solution.
% 
%    See also BETACDF, BETAFIT, BETALIKE, BETAPDF, BETARND, BETASTAT, ICDF.
%
%    Documentation for betainv
%       doc betainv
%
%    Other uses of betainv
%
%       distributed/betainv    rsymbdiff/betainv    splanar/betainv
%