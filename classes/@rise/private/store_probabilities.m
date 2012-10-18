function obj=store_probabilities(obj)
initializing=isempty(obj.state_names);
nstates=sum(max(obj.Regimes,[],1));
if initializing
    obj.state_names=cell(1,nstates);
end
nobs=numel(obj.varobs(1).value);
filt=nan(nstates,nobs+1);
updated=nan(nstates,nobs);
smooth=nan(nstates,nobs);
iter=0;
for ic=1:size(obj.markov_chains,2)
    chain=obj.markov_chains(ic).name;
    for istate=1:size(obj.markov_chains(ic).states,2)
        iter=iter+1;
        regimes=obj.markov_chains(ic).states(istate).regime;
        if initializing
            obj.state_names{iter}=[chain,'_',int2str(istate)];
        end
        if isfield(obj.Filters,'filtered_probabilities')
            filt(iter,:)  =sum(obj.Filters.filtered_probabilities(regimes,:),1);
        end
        if isfield(obj.Filters,'updated_probabilities')
            updated(iter,:)  =sum(obj.Filters.updated_probabilities(regimes,:),1);
        end
        if isfield(obj.Filters,'smoothed_probabilities')
            smooth(iter,:)  =sum(obj.Filters.smoothed_probabilities(regimes,:),1);
        end
    end
end
if obj.options.kf_filtering_level
    obj.Filters.filtered_probabilities=filt;
    if obj.options.kf_filtering_level>1
        obj.Filters.updated_probabilities=updated;
        if obj.options.kf_filtering_level>2
            obj.Filters.smoothed_probabilities=smooth;
        end
    end
end

