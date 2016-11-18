function theStruct=load_priors(theStruct,priors,new_dirichlet)

% for efficiency, this should be done at estimation time?...
if ~isempty(priors)
    % load the distributions
    tmp={priors.prior_distrib};
    
    if ~isempty(tmp)
        
        effective_distributions=unique(tmp);
        
        distr_locs=cell(1,numel(effective_distributions));
        
        for ii=1:numel(effective_distributions)
            
            distr_locs{ii}=find(strcmp(effective_distributions{ii},tmp));
            
            % get the handle on the distributions but not for the
            % dirichlet: they need to be processed separately.
            if ~strcmp(effective_distributions{ii},'dirichlet')
                
                lndens=distributions.(effective_distributions{ii})();
                
                effective_distributions{ii}=lndens;
                
            end
            
        end
        
        theStruct.estim_hyperparams=[[priors.a]',[priors.b]'];
        
        theStruct.estim_distributions=effective_distributions;
        
        theStruct.estim_distrib_locations=distr_locs;
        
    end
    
    theStruct.estim_dirichlet=new_dirichlet;
    
end

end

