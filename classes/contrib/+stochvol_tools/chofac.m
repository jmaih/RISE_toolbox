function CF = chofac(N,chovec);

% CF = chofac(N,chovec);
% This transforms a vector of cholesky coefficients, chovec, into a lower
% triangular cholesky matrix, CF, of size N.

CF = eye(N,N); 
i = 1; 
for j = 2:N,
   k = 1;
   while k < j;
      CF(j,k) = chovec(i,1);
      i = i+1; 
      k = k+1;
   end
end
  