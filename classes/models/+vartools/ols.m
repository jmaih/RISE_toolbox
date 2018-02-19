function out = ols(kdata)

nv=kdata.nvars*kdata.ng;

Y=kdata.estim_.Y;

X=kdata.estim_.X;

if isempty(kdata.estim_.linres)

B=Y/X; % <- Y*X'*inv(X*X')

else
    
    B=use_linear_restrictions();
    
end

Resids=vartools.residuals(kdata,B);
%
T=kdata.estim_.T;

K=kdata.estim_.K;

Sigma=(1/(T-K))*(Resids*Resids.');
%
out=struct('B',B,'Sigma',Sigma);%,'Resids',Resids

    function B=use_linear_restrictions()
        
        Z=kron(X.',speye(nv));
        
        y=Y(:)-Z*kdata.estim_.linres.d;
        
        Z=Z*kdata.estim_.linres.K;
        
        a2tilde=Z\y;
        
        B=reshape(kdata.estim_.linres.a2tilde_to_a(a2tilde),nv,[]);

    end

end
