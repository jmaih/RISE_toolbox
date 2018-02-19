function theStruct=load_priors(theStruct,priors,new_dirichlet,prilocs)

if nargin<4
    
    prilocs=[];
    
end

% for efficiency, this should be done at estimation time?...
if ~isempty(priors)
    % load the distributions
    tmp={priors.prior_distrib};
    
    if ~isempty(tmp)
        
        effective_distributions=unique(tmp);
        
        ned=numel(effective_distributions);
        
        distr_locs=cell(1,ned);
        
        if ~isempty(prilocs)
            
            prilocs=prilocs(:).';
            
        end
        
        for ii=1:ned
            
            distr_locs{ii}=find(strcmp(effective_distributions{ii},tmp));
            
            distr_locs{ii}=distr_locs{ii}(:).';
            
            if isempty(prilocs)
                
                distr_locs{ii}=distr_locs{ii}+distr_locs{ii}*1i;
                
            else
                
                distr_locs{ii}=distr_locs{ii}+prilocs(distr_locs{ii})*1i;
                
            end
            
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

