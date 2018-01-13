function [Va,V,sig] = set_prior_variance(Yraw,SIGu,kdata,prior_hyperparams)

V=[];

sig=[];

switch lower(prior_hyperparams.type)
    
    case 'minnesota'
        
        Va = set_minnesota_variance();
        
    case {'normal-wishart','nw'}
        
        Va = set_normal_wishart_variance();
        
    case {'indep-normal-wishart','inw'}
        
        Va = set_independent_normal_wishart_variance();
        
    case {'sims-zha','sz'}
        
        Va = [];
        
        sig=prior_sigma();
        
    otherwise
        
        error(['unknown prior type ',prior_hyperparams.type])
end

    function Va = set_independent_normal_wishart_variance()
        
        if isscalar(prior_hyperparams.independent_normal_wishart_eta)
            
            Va=diag(prior_hyperparams.independent_normal_wishart_eta*ones(1,kdata.K*kdata.nvars));
            
        elseif isequal(size(prior_hyperparams.normal_wishart_eta),kdata.K*kdata.nvars*ones(1,2))
            
            Va=prior_hyperparams.independent_normal_wishart_eta;
            
        else
            
            error('wrong format for the independent normal wishart variance')
            
        end
        
    end

    function Va = set_normal_wishart_variance()
        
        if isscalar(prior_hyperparams.normal_wishart_eta)
            
            V=diag(prior_hyperparams.normal_wishart_eta*ones(1,kdata.K));
            
        elseif isequal(size(prior_hyperparams.normal_wishart_eta),kdata.K*ones(1,2))
            
            V=prior_hyperparams.normal_wishart_eta;
            
        else
            
            error('wrong format for the normal wishart within equation variance')
            
        end
        
        Va=kron(V,SIGu); % 5.2.13
        
    end

    function [Va,V] = set_minnesota_variance()
        
        V=[];
        
        L1=prior_hyperparams.L1;
        
        L2=prior_hyperparams.L2;
        
        L3=prior_hyperparams.L3;
        
        L4=prior_hyperparams.L4;
        
        Sigma=prior_sigma();
        
        Va = zeros(kdata.nvars,kdata.K);
        
        Va(:,1:kdata.nx) = (Sigma(:)*ones(1,kdata.nx))*(L1*L4)^2;
        
        for ll = 1:kdata.nlags
            
            for jrow = 1:kdata.nvars
                
                for kcol = 1:kdata.nvars
                    
                    Va(jrow,(ll-1)*kdata.nvars + kcol + kdata.nx) = ...
                        lag_variance(ll,jrow,kcol);
                    
                end
                
            end
            
        end
        
        Va = diag(Va(:));
        
        function v=lag_variance(ilag,jrow,kcol)
            
            v=(Sigma(jrow)/Sigma(kcol))*...
                (L1/(ilag^L3))^2;
            
            if jrow~=kcol
                
                v=v*L2^2;
                
            end
            
        end
        
    end
        
    function s=prior_sigma()
        
        s=zeros(1,kdata.nvars);
        
        % note Yraw could be a panel... in which case it will have nt
        % columns and np pages. Hence nt_np is the product of both
        [~,nt_np]=size(Yraw);
        
        % add a constant to the regression
        %----------------------------------
        x=ones(1,nt_np);
            
        kdatai=struct('nlags',kdata.nlags,'ng',1,...
            'nvars',1,'nx',1,'linres',[]);
        
        for iii = 1:kdata.nvars
            % expand for (homogenous) panel...
            kdatai = rfvar.embed(kdatai,Yraw(iii,:),x);
            
            mlei=vartools.ols(kdatai);
            
            s(iii) = mlei.Sigma;
            
        end
        
    end

end