function p=param_creator(nvars,lags,constant,nd,additional_params,prefix,is_svar)

% constant is the number of constants, which could be more than 1 in the
% case of a panel!!!

lags=lags(1):lags(end);

nlags=numel(lags);

% np=(nlags*nvars+nx)*nvars;

pnames=cell(nvars,nlags*nvars+nd+constant+is_svar+(1-is_svar)*nvars);

shadow=nan(nvars,nlags*nvars+nd+constant+is_svar+(1-is_svar)*nvars);

lb=-inf(size(shadow));

ub=-lb;

start=zeros(size(shadow));

for eqtn=1:nvars
    
    iter=0;
    
    % deterministic
    for ix=1:nd+constant
        
        iter=iter+1;
        
        pnames{eqtn,iter}=sprintf('c_%0.0f_%0.0f',eqtn,ix);
        
    end
    
    % lags
    for ilag=lags
        
        for iv=1:nvars
            
            iter=iter+1;
            
            pnames{eqtn,iter}=sprintf('%s%0.0f_%0.0f_%0.0f',prefix,ilag,eqtn,iv);
            
            if is_svar && ilag==0
                % set the diagonal of the shadow to 1: the parameters won't
                % be estimated. Instead the standard deviations on the
                % right-hand side will be.
                if iv==eqtn
                    
                    shadow(eqtn,iter)=1;
                    
                end
                
            end
            
        end
        
    end
    
    if is_svar
        % standard deviations for each equation
        pnames{eqtn,iter+1}=sprintf('s_%0.0f_%0.0f',eqtn,eqtn);
        
        lb(eqtn,iter+1)=sqrt(eps);
        
        start(eqtn,iter+1)=0.5;
        
    else
        % covariance matrix
        for iv=1:nvars
            % take the vech
            if iv>eqtn
                
                continue
                
            end
            
            iter=iter+1;
            
            pnames{eqtn,iter}=sprintf('s_%0.0f_%0.0f',eqtn,iv);
            
            if eqtn==iv
                
                lb(eqtn,iter)=sqrt(eps);
                
                start(eqtn,iter)=0.5;
                
            end
            
        end
        
    end
    
end

pnames=pnames(:);

lb=lb(:);

ub=ub(:);

shadow=shadow(:);

start=start(:);

% remove empty entries
bad=cellfun(@isempty,pnames);

pnames(bad)=[];

lb(bad)=[];

ub(bad)=[];

shadow(bad)=[];

start(bad)=[];

% adding additional parameters (transition probabilities and others)
%-------------------------------------------------------------------
nadd=numel(additional_params);

if nadd
    
    pnames=[pnames(:);additional_params(:)];
    
    lb=[lb(:);nan(nadd,1)];
    
    ub=[ub(:);nan(nadd,1)];
    
    shadow=[shadow(:);nan(nadd,1)];
    
    start=[start(:);nan(nadd,1)];
    
end

p=struct('pnames',{pnames},'shadow',shadow,'lb',lb,'ub',ub,'start',start,...
    'var_terms',nvars*(nlags*nvars+nd+constant),...
    's_terms',nvars*is_svar+(1-is_svar)*sum(1:nvars),...
    'tp_terms',sum(parser.is_transition_probability(pnames)));

end