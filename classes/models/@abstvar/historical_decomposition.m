function [hd,retcode]=historical_decomposition(self,params,Rfunc,shock_names)
% varargin={shock_names}

n=nargin;

set_defaults()

np=numel(params);

retcode=nan(1,np);

for ip=1:np
    
    if ip==1
        
        [hd,retcode(ip),info,contrib_pos]=do_one_parameter(params(:,ip)); %#ok<ASGLU>
        
        hd=hd(:,:,:,ones(1,np));
        
    else
        
        [hd(:,:,:,ip),retcode(ip)]=do_one_parameter(params(:,ip));
        
    end
    
end

tmp=struct();

start_date=self.estim_.date_range(1)+self.nlags;

contrib_names=cell(1,size(hd,2));

contrib_names(contrib_pos.shocks)=abstvar.create_variable_names(self.nvars,'shock',shock_names);

contrib_names{contrib_pos.y0}='y0';

contrib_names(contrib_pos.det_vars)=abstvar.create_variable_names(numel(contrib_pos.det_vars),'det',[]);

for iv=1:self.nvars
    
    tmp.(self.endogenous{iv})=ts(start_date,permute(hd(:,:,iv,:),[1,2,4,3]),contrib_names);
    
end

hd=tmp;

    function [hd,retcode,info,contrib_names]=do_one_parameter(param)
        
        Resids=vartools.residuals(self,param.B);
        
        [R,retcode]=Rfunc(param);
        
        [hd,info,contrib_names]=vartools.historical_decomposition(param.B,...
            R,self.nx*self.ng,self.estim_.X,Resids);
        
    end

    function set_defaults()
        
        if n<4
            
            shock_names=[];
            
            if n<3
                
                Rfunc=double.empty;
                
                if n<2
                    
                    params=[];
                    
                end
                
            end
            
        end
        
        params=solve(self,params);
        
        if isempty(Rfunc),Rfunc=identification(self,'choleski'); end
        
    end

end

