function filtering=store_probabilities(obj,filtering)
% INTERNAL FUNCTION
%

nstates=numel(obj.markov_chains.state_names);

nobs=obj.data.finish;

filt=zeros(nstates,nobs+1);

updated=zeros(nstates,nobs);

smooth=zeros(nstates,nobs);

iter=0;

regimes=cell2mat(obj.markov_chains.regimes(2:end,2:end));

for ic=1:obj.markov_chains.chains_number
    
    max_states=max(regimes(:,ic));
    
    for istate=1:max_states
        
        iter=iter+1;
        
        this_regimes=find(regimes(:,ic)==istate);
        
        for ireg=1:numel(this_regimes)
            
            if isfield(filtering,'filtered_regime_probabilities')
                
                filt(iter,:)=filt(iter,:)+filtering.filtered_regime_probabilities(this_regimes(ireg),:);
            
            end
            
            if isfield(filtering,'updated_regime_probabilities')
                
                updated(iter,:)=updated(iter,:)+filtering.updated_regime_probabilities(this_regimes(ireg),:);
            
            end
            
            if isfield(filtering,'smoothed_regime_probabilities')
                
                smooth(iter,:)=smooth(iter,:)+filtering.smoothed_regime_probabilities(this_regimes(ireg),:);
            
            end
            
        end
        
    end
    
end

if obj.options.kf_filtering_level
    
    filtering.filtered_state_probabilities=filt;
    
    if obj.options.kf_filtering_level>1
        
        filtering.updated_state_probabilities=updated;
        
        if obj.options.kf_filtering_level>2
            
            filtering.smoothed_state_probabilities=smooth;
            
        end
        
    end
    
end

