function RepsRun=bootstrap(self,boot)

if nargin < 2
    
    boot=[];
    
end

if isempty(boot),boot=1000; end

if ~isempty(self.prior)
    
    warning('No bootstrapping for models estimated using Bayesian techniques')
    
    if (self.is_switching||self.optimize)
        
        error('no bootstrapping for models estimated using bayesian techniques or optimzed numerically')
        
    else
        
        warning('switching to analytical sampling algorithms for constant-parameter BVARs')
        
        RepsRun=self.sampler(boot);
        
    end
    
end

sol=solve(self);

B=sol.B;

Resids=self.Y(:,:)-B*self.X(:,:);

nv=self.nvars;

nt=self.T;

ng=self.ng;

Y=self.Y;

X=self.X;

T=self.T;

% K=kdata.K;

ndet=self.nx;

nvars=size(Y,1);

population=1:T;
%
fprintf('Starting bootstrapping %s\n', datestr(now));
%
for ii =1:boot
    
    SPoints = randsample(population,T,true);
    
    if ii==1
        
        RepsRun = cvarreps(Resids(:,SPoints));
        
        RepsRun=RepsRun(:,ones(boot,1));
        
    else
        
        RepsRun(:,ii) = cvarreps(Resids(:,SPoints));
        
    end
    
end

fprintf('bootstrapping reps finished, %s\n', datestr(now));

    function RepsRun = cvarreps(Resids)
        
        Yi=Y;
        
        Xi=X;
        
        for jcol=1:T
            
            Yi(:,jcol)=B*Xi(:,jcol)+Resids(:,jcol);
            
            if jcol<T
                
                Xi(ndet+1:end,jcol+1)=[Yi(:,jcol)
                    Xi(ndet+1:end-nvars,jcol)];
                
            end
            
        end
        
        kdatai=update_base_data(self,Xi,Yi);
        
        RepsRun = vartools.ols(kdatai);
        
        RepsRun=[RepsRun.B(:);vech(RepsRun.Sigma)];
        
        function k=update_base_data(k,X,Y)
            
            k.X=X;
            
            k.Y=Y;
            
        end
        
    end

end

