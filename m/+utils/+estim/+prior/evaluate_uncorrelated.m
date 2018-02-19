function [lnprior,retcode]=evaluate_uncorrelated(epdata,param,lnprior)

retcode=0;

if isempty(epdata)
    
    return
    
end

for kk=1:numel(epdata.estim_distributions)
    
    if ischar(epdata.estim_distributions{kk})
        
        % that is the dirichlet,skip it
        continue
        
    end
    
    loc=epdata.estim_distrib_locations{kk};
    
    % location in the structure itself
    loc1=real(loc);
    
    % location in the grand vector of all estimated parameters
    loc2=imag(loc);
    
    a=epdata.estim_hyperparams(loc1,1);
    
    b=epdata.estim_hyperparams(loc1,2);
    
    ld=epdata.estim_distributions{kk}(param(loc2),a,b);
        
    lnprior=lnprior+sum(ld);
    
    if (isnan(lnprior)||isinf(lnprior)||~isreal(lnprior))
        
        retcode=307;
        
        break
        
    end
    
end

end
