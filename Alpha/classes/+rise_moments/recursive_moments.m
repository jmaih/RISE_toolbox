function [m1,V1]=recursive_moments(m0,V0,Xn,n)
m1=1/n*(Xn+(n-1)*m0);
if nargout>1
    V1=(n-1)/n*(V0+m0*m0')+1/n*(Xn*Xn')-m1*m1';
end

%{
npar=11;
ndraws=3000;
Params=rand(npar,ndraws);
m=0;
V=0;
for ii=1:ndraws
    [m,V]=recursive_moments(m,V,Params(:,ii),ii);
end
max(max(abs(V-cov(Params',1))))
%}