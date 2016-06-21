function SIG=initial_covariance(lb,ub,k)

if nargin<3
    
    k=5;
    
end

ub(ub==inf)=10;

lb(lb==-inf)=-10;

sd=(ub-lb)/(2*k);

SIG=diag(sd.^2); % SIG=1e-4*eye(d);%SIG=covar;

end