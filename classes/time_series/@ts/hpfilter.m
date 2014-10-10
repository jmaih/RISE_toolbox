function [hptrend,hpcycle] = hpfilter(db,lambda)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

% HP filters a collection of time series.
% 
% INPUTS 
%   db                       [double]   T*n ts object
%   lambda                   [double]   scalar, lambda parameter.
% 
% OUTPUTS 
%   hptrend                  [double]   T*n ts object, trend component of db.
%   hpcycle                  [double]   T*n ts object, cycle component of db.  
%               

y=db.data;
T = db.NumberOfObservations;

if nargin<2
    lambda = 1600;
end

A=repmat([lambda, -4*lambda, (1 + 6*lambda), -4*lambda, lambda], T, 1);
A = spdiags(A, -2:2, T,T);
A(1,1) = 1+lambda; 
A(1,2) = -2*lambda; 
A(2,2) = 1+5*lambda; 
A(T-1,T-1) = 1+5*lambda;
A(T-1,T) = -2*lambda; 
A(T,T-1) = -2*lambda;
A(2,1) = -2*lambda; 
A(T,T) = 1+lambda;

hptrend = ts(db.start,A\y,db.varnames);
hpcycle = db-hptrend;
hpcycle.varnames=db.varnames;

end


