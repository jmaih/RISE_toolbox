function [f,the_regimes,endog,data,tex]=load_filters(m)

n=numel(m);

if n>1
    
    f=cell(size(m));
    the_regimes=f;
    endog=f;
    data=f;
    tex=f;
    
    for ii=1:n
        
        [f{ii},the_regimes{ii},endog{ii},data{ii},tex{ii}]=load_filters(m(ii));
        
    end
    
    return
    
end

if isa(m,'abstvar')
    
    [~,~,~,f]=filter(m);
    
    the_regimes=generic.describe_regimes(m.markov_chain_info);
    
    data=m.estim_.data;
    
    endog=m.endogenous;
    
    tex=[];
    
else
    
    f=filter(m);
    
    the_regimes=generic.describe_regimes(m.markov_chains);
    
    endog=m.observables.name;
    
    data=m.options.data;
    
    if isa(data,'ts')
        
        data=pages2struct(data);
        
    end
    
    tex=get(m,'tex');
    
end

end