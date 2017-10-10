function [Y,X,Yc,Xc]=sims_zha_dummies(prior,nvars,nx,nlags,sig,ybar)

B=vartools.bvar_coef_prior(prior.coefprior,nvars,nlags,nx);

B1_ii=diag(B(:,nx+(1:nvars)));

% L1 : Overall tightness
tau=1/prior.L1;

% L2 : Cross-variable specific variance parameter

% L3 : Speed at which lags greater than 1 converge to zero
d=1/prior.L3;

% L4 : tightness on deterministic/exogenous terms

% covariance dummies
omega=prior.L5;

if ~(isscalar(omega) && isnumeric(omega) && ...
        isfinite(omega) && floor(omega)==ceil(omega))
    
    error('L5 must be a finite integer')
    
end

% co-persistence
lambda=prior.L6;

% Own-persistence
mu=prior.L7;

dtypes={@lag_dummies,@covariance_dummies,...
    @copersistence_dummies,@own_persistence_dummies};

ntypes=numel(dtypes);

Yc=cell(1,ntypes);

Xc=cell(1,ntypes);

for itype=1:ntypes
    
    [Yc{itype},Xc{itype}]=dtypes{itype}();
    
end

Y=cell2mat(Yc);

X=cell2mat(Xc);

    function [Y,X]=lag_dummies()
        
        Y=cell(1,nlags);
        
        X=cell(1,nlags);
        
        Y0=zeros(nvars,nvars);
        
        X0=zeros(nx+nvars*nlags,nvars);
        
        tmp=sig*tau;
        
        for ilag=1:nlags
            
            Yi=Y0; Xi=X0;
            
            if ilag==1
                
                Yi=diag(tmp.*B1_ii);
                
            end
            
            Xi(nx+(ilag-1)*nvars+(1:nvars),:)=diag(tmp*ilag^d);
            
            Y{ilag}=Yi;
            
            X{ilag}=Xi;
            
        end
        
        Y=cell2mat(Y);
        
        X=cell2mat(X);

    end

    function [Y,X]=covariance_dummies()
        
        Y0=diag(sig);
        
        X0=zeros(nx+nvars*nlags,nvars);
        
        Y=cell(1,nlags);
        
        X=cell(1,nlags);
        
        for ii=1:omega
            
            Y{ii}=Y0;
            
            X{ii}=X0;
            
        end
        
        Y=cell2mat(Y);
        
        X=cell2mat(X);
        
    end

    function [Y,X]=copersistence_dummies()
        
        lambda_negative=lambda<0;
        
        if lambda_negative
            
            lambda=-lambda;
            
        end
        
        ly=lambda*ybar(:);
        
        Y=ly;
        
        ly=ly(:,ones(nlags,1));
        
        X=[(1-lambda_negative)*lambda*ones(nx,1)
            ly(:)];
        
    end

    function [Y,X]=own_persistence_dummies()
        
        tmp=diag(mu*ybar(:));
        
        Y=tmp;
        
        X=zeros(nx+nvars*nlags,nvars);
                
        for ilag=1:nlags
            
            X(nx+(ilag-1)*nvars+(1:nvars),:)=tmp;
            
        end
        
    end

end