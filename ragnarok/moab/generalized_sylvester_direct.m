function [Tzz,retcode]=generalized_sylvester_direct(A,B,C,fv_lamba_bf,W)

need_uncompress = nargin==5;

[nt,nzz]=size(A{1});

h=numel(A);

G=zeros(nt*nzz*h);

offsetr=0;

a=zeros(nt*nzz*h,1);

for rt=1:h
    
    r_range=offsetr+(1:nt*nzz);
    
    a(r_range)=A{rt}(:);
    
    offsetc=0;
    
    for rtp1=1:h
        
        c_range=offsetc+(1:nt*nzz);
        
        G(r_range,c_range)=kron(C{rt,rtp1}.',fv_lamba_bf{rt,rtp1});
        
        if rt==rtp1
            
            G(r_range,c_range)=G(r_range,c_range)+kron(speye(nzz),B{rt});
            
        end
        
        offsetc=c_range(end);
        
    end
    
    offsetr=r_range(end);
    
end

G=sparse(G);

Tzz0=-G\a;

retcode=1-all(isfinite(Tzz0));

Tzz0=reshape(Tzz0,nt,nzz,h);

Tzz=cell(1,h);

for rt=1:h
    
    if need_uncompress
        
        Tzz{rt}=sparse(Tzz0(:,:,rt)*W);
        
    else
        
        Tzz{rt}=sparse(Tzz0(:,:,rt));
        
    end
    
end

end
