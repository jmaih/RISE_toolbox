function [C,R]=autocov(self,params,max_periods)

if nargin<3
    
    max_periods=[];
    
    if nargin<2
        
        params=[];
        
    end
    
end

if isempty(max_periods),max_periods=5; end

params=solve(self,params);

Rfunc=identification(self,'choleski');

np=numel(params);

for ip=1:np
    
    if ip==1
        
        [C,R,info]=do_one_parameter(params(:,ip)); %#ok<ASGLU>
        
        C=C(:,:,:,ones(1,np));
        
        R=R(:,:,:,ones(1,np));
        
    else
        
        [C(:,:,:,ip),R(:,:,:,ip)]=do_one_parameter(params(:,ip));
        
    end
    
end

    function [C,R,info]=do_one_parameter(param)
        
        [C,R,info]=vartools.autocorr(param.B,...
            Rfunc(param),...
            self.nx*self.ng,max_periods);
        
    end

end

