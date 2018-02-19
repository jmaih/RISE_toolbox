function [f,the_regimes,endog,data,tex]=load_filters(m)

if isa(m,'abstvar')
    
    [~,~,~,f]=filter(m);
    
    the_regimes=generic.describe_regimes(m.markov_chain_info);
    
    data=m.estim_.data;
    
    endog=m.endogenous;
    
    tex=[];
    
else
    
    f=m.filtering;
    
    the_regimes=generic.describe_regimes(m.markov_chains);
    
    endog=m.observables.name;
    
    data=m.options.data;
    
    if isa(data,'ts')
        
        data=pages2struct(data);
        
    end
    
    tex=get(m,'tex');
    
end

end