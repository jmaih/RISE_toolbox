function [hptrend,hpcycle]=hpfilter(y,lambda)

if nargin<2
    
    lambda = 1600;
    
end

T = size(y,1);

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

hptrend = A\y;

hpcycle = y-hptrend;

end