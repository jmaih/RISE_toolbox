function out = ols(kdata)

nv=kdata.nvars*kdata.ng;

Y=kdata.Y;

X=kdata.X;

if isempty(kdata.linres)

B=Y/X; % <- Y*X'*inv(X*X')

else
    
    B=use_linear_restrictions();
    
end

Resids=vartools.residuals(kdata,B);
%
T=kdata.T;

K=kdata.K;

Sigma=(1/(T-K))*(Resids*Resids.');
%
out=struct('B',B,'Sigma',Sigma);%,'Resids',Resids

    function B=use_linear_restrictions()
        
        Z=kron(X.',speye(nv));
        
        y=Y(:)-Z*kdata.linres.d;
        
        Z=Z*kdata.linres.K;
        
        a2tilde=Z\y;
        
        B=reshape(kdata.linres.a2tilde_to_a(a2tilde),nv,[]);

    end

end
