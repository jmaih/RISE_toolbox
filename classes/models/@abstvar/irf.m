function [myirfs]=irf(self,shock_names,irf_periods,params,Rfunc)

% N.B: Only non-failed irfs are returned!!!

n=nargin;

set_defaults()

[myirfs,info]=vartools.irf(self,irf_periods,params,Rfunc); %#ok<ASGLU>

tmp=struct();

shock_names=abstvar.create_variable_names(self.nvars,'shock',shock_names);

for ishock=1:self.nvars
    
    for iv=1:self.nvars
        
        data=permute(myirfs(iv,:,ishock,:,:),[2,5,4,1,3]);
        
        data=squeeze(data);
        
        if ishock==1 && iv==1
            
            if self.nregs==1
                
                proto=ts(1,data);
                
            else
                
                regimeNames=abstvar.create_variable_names(self.nregs,'regime');
                
                proto=ts(1,data,regimeNames);
                
            end
            
        end
        
        tmp.(shock_names{ishock}).(self.endogenous{iv})=reset_data(proto,data);
        
    end
    
end

myirfs=tmp;

    function set_defaults()
        
        if n < 5
            
            Rfunc=[];
            
            if n < 4
                
                params=[];
                
                if n< 3
                    
                    irf_periods=[];
                    
                    if n<2
                        
                        shock_names=[];
                        
                    end
                    
                end
                
            end
            
        end
        
        if isempty(Rfunc),Rfunc=identification(self,'choleski'); end
        
        params=solve(self,params);
        
        if isempty(irf_periods),irf_periods=40; end
        
    end

end

