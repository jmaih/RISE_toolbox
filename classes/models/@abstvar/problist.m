function [probnames,nprobs]=problist(cn,nstates)

nprobs=nstates^2-nstates;

probnames=cell(1,nprobs);

iter=0;

for ii=1:nstates
    
    for jj=1:nstates
        
        pname=sprintf('%s_tp_%0.0f_%0.0f',cn,ii,jj);
        
        if ii==jj,continue,end
        
        iter=iter+1;
        
        probnames{iter}=pname;
        
    end
    
end

end