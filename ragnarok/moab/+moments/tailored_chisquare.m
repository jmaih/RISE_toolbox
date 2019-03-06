function  m=tailored_chisquare(order,k)
% Inspired by a question on the dynare forum at
% http://www.dynare.org/phpBB3/viewtopic.php?f=1&t=3206 

% moments of x=(e-k)/sqrt(2*k)

 m0=moments.chisquare(order,k);
 
 m=zeros(order,1);
 
 m0=[1,m0(:).']; % add 1 at the beginning so that we have E(e^0)
 
 s2k=sqrt(2*k);
 
 s2ki=1;
 
 pt=moments.pascal_triangle(order);
 
 for io=1:order
     
     pt_i=pt{io+1};
     
     mi=0;
     
     n=numel(pt_i);
     
     for icol=1:n % == order+1
         
         c=pt_i(icol);
         
         mi=mi+c*m0(icol)*(-k)^(n-icol);
         
     end
     
     s2ki=s2ki*s2k;
     
     m(io)=mi/s2ki;
     
 end

end