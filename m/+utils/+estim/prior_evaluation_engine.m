function [lnprior,retcode]=prior_evaluation_engine(epdata,param,lnprior)

if nargin<3
    
    lnprior=0;
    
end

% do the uncorrelated guys
%--------------------------
[lnprior,retcode]=utils.estim.prior.evaluate_uncorrelated(epdata,param,...
    lnprior);

if isempty(epdata)
    
    return
    
end

% do the dirichlet guys
%-----------------------
if ~retcode && ~isempty(epdata.estim_dirichlet)
    
    ndirich=numel(epdata.estim_dirichlet);
    
    for id=1:ndirich
        
        loc=epdata.estim_dirichlet(id).location;
        
        lnprior=lnprior+epdata.estim_dirichlet(id).lpdfn(param(loc));
        
        if (isnan(lnprior)||isinf(lnprior)||~isreal(lnprior))
            
            retcode=307;
            
            break
            
        end
        
    end
    
end


end