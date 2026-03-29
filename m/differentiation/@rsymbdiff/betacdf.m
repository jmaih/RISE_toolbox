% BETACDF Beta cumulative distribution function.
%    P = BETACDF(X,A,B) returns the beta cumulative distribution
%    function with parameters A and B at the values in X.
% 
%    The size of P is the common size of the input arguments. A scalar input  
%    functions as a constant matrix of the same size as the other inputs.    
% 
%    BETAINC does the computational work.
% 
%    P = BETACDF(X,A,B,'upper') returns the upper tail probability of the beta 
%    distribution function with parameters A and B at the values in X.
% 
%    See also BETAFIT, BETAINV, BETALIKE, BETAPDF, BETARND, BETASTAT, CDF,
%             BETAINC.
%
%    Documentation for betacdf
%       doc betacdf
%
%    Other uses of betacdf
%
%       distributed/betacdf    rsymbdiff/betacdf    splanar/betacdf
%