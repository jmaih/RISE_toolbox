function M=estim2states(params,theMap,nparams,nregs)

nsol=size(params,2);

M=zeros(nparams,nregs,nsol);

for isol=1:nsol

    M(:,:,isol)=do_one(M(:,:,isol),params(:,isol));

end

    function M=do_one(M,p)
        
        M(theMap(:,2))=p(theMap(:,1));
        
    end

end
