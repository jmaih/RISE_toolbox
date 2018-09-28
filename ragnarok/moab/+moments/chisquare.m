function m=chisquare(order,k)
% recursive computation of the moments of the chi-square distribution with
% k degrees of freedom
% - order : maximum order for which one wants the moments to be computed
% for
% m : vector of moments from 1 to order
m=zeros(order,1);

m(1)=k;

for o=2:order
    
    m(o)=2*(0.5*k+(o-1))*m(o-1);
    
end

end
