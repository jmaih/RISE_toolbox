function flag=is_has_solution(obj,pnames,pvals)

p=struct();

for ii=1:numel(pnames)
    
    p.(pnames{ii})=pvals(ii);
    
end

[~,retcode]=solve(obj,'parameters',p);

flag=retcode==0;   

end

