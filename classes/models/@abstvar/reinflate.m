function p=reinflate(p,chain_name,nstates)

if ischar(p),p=cellstr(p); end

if strcmp(chain_name,'const')
    
    return
    
end

p=p(:);

p=p(:,ones(1,nstates));

for irow=1:size(p,1)
    
    for istate=1:nstates
        
        p{irow,istate}=sprintf('%s_%s_%0.0f',p{irow,istate},chain_name,istate);
        
    end
    
end

p=p(:).';

end
