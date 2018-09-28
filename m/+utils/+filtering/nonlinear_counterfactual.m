function [C,CD]=nonlinear_counterfactual(y0,fsim,shocks,nsim,varargin)

do_distrib=nargout>1;

if ~do_distrib
    
    C=0;
    
end

for jsim=1:nsim
    
    Ci=fsim(y0,shocks,varargin{:});
    
    if do_distrib
        
        if jsim==1
            
            CD=Ci(:,:,ones(1,nsim));
            
        else
            
            CD(:,:,jsim)=Ci;
            
        end
        
    else
        
        C=C+Ci;
        
    end
    
end

if do_distrib
    
    C=mean(CD,3);
    
else
    
    C=1/nsim*C;
    
end

end