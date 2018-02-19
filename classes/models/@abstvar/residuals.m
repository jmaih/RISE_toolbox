function [r,f]=residuals(self,params,shock_names,raw)

n=nargin;

set_defaults()

np=numel(params);

ng=self.ng;

if ng==1
    
    countryNames={'fake'};
    
else
    
    countryNames=self.members;
    
end

tmpres=struct();

tmpfit=struct();

for ip=1:np
    
    [rx,fx]=do_one_parameter(params(:,ip));
    
    if ip==1
        
        nvars=size(rx,1)/ng;
        
    end
    
    for ig=1:ng
        
        batch=ig:ng:nvars*ng;
        
        store_one_parameter_country(rx(batch,:,:),fx(batch,:,:))
        
    end
    
end

start_date=serial2date(self.estim_.date_range(1)+self.nlags);

if ng==1
    
    tmpres=tmpres.(countryNames{1});
    
    tmpfit=tmpfit.(countryNames{1});
    
end

if raw
    
    r=struct('residuals',tmpres,'start_date',start_date);
    
    f=struct('fitted',tmpfit,'start_date',start_date);
    
    return
    
end

r=set_data_to_time_series(self,tmpres,shock_names,start_date);

f=set_data_to_time_series(self,tmpfit,self.endogenous,start_date);

    function store_one_parameter_country(rx,fx)
        
        if ip==1
            
            tmpres.(countryNames{ig})=rx(:,:,:,ones(1,np));
        
            tmpfit.(countryNames{ig})=fx(:,:,:,ones(1,np));
        
        else
            
            tmpres.(countryNames{ig})(:,:,:,ip)=rx;
        
            tmpfit.(countryNames{ig})(:,:,:,ip)=fx;
        
        end
        
    end

    function [Resids,Fits]=do_one_parameter(param)
        
        [Resids,Fits]=vartools.residuals(self,param.B);
        
    end

    function set_defaults()
        
        if n<4
            
            raw=[];
            
            if n<3
                
                shock_names=[];
                
                if n<2
                    
                    params=[];
                    
                end
                
            end
            
        end
        
        if isempty(shock_names)
            
            shock_names=abstvar.create_variable_names(self.nvars,'shock',shock_names);
            
        end
        
        if isempty(raw)
            
            raw=false;
            
        end
        
        params=solve(self,params);
        
    end

end

