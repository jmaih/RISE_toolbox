function flag=price_puzzle(obj,pnames,pvals)

p=struct();

for ii=1:numel(pnames)
    
    p.(pnames{ii})=pvals(ii);
    
end

[obj,retcode]=solve(obj,'parameters',p);

if retcode
    
    flag=false;
    
else
    
    [state_list]=state_var_list(obj,1);
    
    % inflation
    vloc=strcmp('PIE',obj.endogenous.name);
    
    % monetary policy shock
    sloc=strcmp('EI',state_list);
    
    response=obj.state_space{1}.Tz{1}(vloc,sloc);
    
    flag=response>0;
    
end

end

