function [sol,M,retcode]=solve(mapping,param)

nregs=size(mapping.regimes,1)-1;

nparams=mapping.nparams;

a=0.01;

np=size(param,2);

retcode=zeros(1,np);

nes=numel(mapping.estimList);

if size(param,1)~=nes
    
    error('wrong size of the parameter input')
    
end

theMap=mapping.theMap;

transProbs=mapping.transProbs;

nprob=numel(transProbs);

for r=1:np
    
    tmpM=estim2states(param(:,r));
    
    [tmpSol,retcode(r)]=mapping.solution(tmpM);
    
    if r==1
        
        M=tmpM(:,:,ones(np,1));
        
        sol=tmpSol(1,ones(np,1));
        
    end
    
    M(:,:,r)=tmpM;
    
    sol(r)=tmpSol;
    
end

    function M=estim2states(param)
        
        for ii=1:nprob
            
            loc=transProbs{ii};
            
            param(loc)=reprobabilize(param(loc));
            
        end
        
        M=zeros(nparams,nregs);
        
        M(theMap(:,2))=param(theMap(:,1));
        
    end

    function x=reprobabilize(x)
        
        x2=x.^2;
        
        x=x2/(a+sum(x2));
        
    end

end