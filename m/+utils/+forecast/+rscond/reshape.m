function shocks=reshape(A,nshocks,debug)
[nr,np]=size(A);
nc=nr/nshocks;
shocks=reshape(A,[nshocks,nc,np]);
if debug
    shocks2=zeros(nshocks,nc,np);
    for ipage=1:np
        shocks2(:,:,ipage)=reshape(A(:,ipage),nshocks,nc);
    end
    shocks3=zeros(nshocks,nc,np);
    offset=0;
    for ipage=1:np
        shocks3(offset+(1:nr))=A(:,ipage);
        offset=offset+nr;
    end
    disp('checking reshape shocks')
    all(all(all(shocks-shocks2==0)))
    all(all(all(shocks-shocks3==0)))
end
end