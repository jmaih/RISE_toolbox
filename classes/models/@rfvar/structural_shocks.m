function [s,retcode]=structural_shocks(self,params,Rfunc,shock_names)

nargs=nargin;

set_defaults()

s=residuals(self,params,shock_names,true);

start_date=s.start_date;

s=s.residuals;

is_panel=self.is_panel;

if is_panel
    
    members=self.members;
    
    ng=numel(members);
    
end

np=numel(params);

retcode=nan(1,np);

nregs=size(params(1).B,3);

nvars=self.nvars;

ng=self.ng;

for ip=1:np
    
    [R,retcode(ip)]=Rfunc(params(ip));
        
    for ireg=1:nregs
        
        iR=panelize(R(:,:,ireg)\eye(nvars*ng));
        
        if is_panel
            
            for g=1:ng
                
                s.(members{g})(:,:,ireg,ip)=iR(:,:,g)*s.(members{g})(:,:,ireg,ip);
                
            end
            
        else
            
            s(:,:,ireg,ip)=iR*s(:,:,ireg,ip);
                
        end
        
    end
    
end

s=set_data_to_time_series(self,s,shock_names,start_date);

    function iR=panelize(iR)
        
        if ng==1
            
            return
            
        end
        
        tmp=iR;
        
        iR=zeros(nvars,nvars,ng);
        
        for gg=1:ng
            
            batch_rows=(gg-1)*nvars+1:gg*nvars;
            
            iR(:,:,gg)=tmp(batch_rows,batch_rows);
            
        end
        
    end

    function set_defaults()
        
        if nargs < 4
            
            shock_names=[];
            
            if nargs < 3
                
                Rfunc=double.empty;
                
                if nargs<2
                    
                    params=[];
                    
                end
                
            end
            
        end
        
        if isempty(Rfunc),Rfunc=identification(self,'choleski'); end
        
        params=solve(self,params);
        
        shock_names=abstvar.create_variable_names(self.nvars,'shock',shock_names);
        
    end

end
