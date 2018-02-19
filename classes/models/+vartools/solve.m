function [sol,retcode]=solve(mapping,M,yt)

nsol=size(M,3);

retcode=zeros(1,nsol);

for r=1:nsol
    
    [tmpSol,retcode(r)]=mapping.solution.full(M(:,:,r),yt);
    
    if r==1
               
        sol=tmpSol(1,ones(nsol,1));
        
    end
    
    sol(r)=tmpSol;
    
end


end
