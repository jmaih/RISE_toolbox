function [Tzz,retcode]=generalized_sylvester_direct2(A,B,C,fv_lamba_bf,W)

need_uncompress = nargin==5;

[nt,nzz]=size(A{1});

h=numel(A);

offsetr=0;

a=zeros(nt*nzz*h,1);

ii=cell(h^2,1);

jj=cell(h^2,1);

ss=cell(h^2,1);

iter=0;

for rt=1:h
    
    r_range=offsetr+(1:nt*nzz);
    
    a(r_range)=-A{rt}(:);
    
    offsetc=0;
    
    for rtp1=1:h
        
        iter=iter+1;
        
        c_range=offsetc+(1:nt*nzz);
        
        Grc=kron(C{rt,rtp1}.',fv_lamba_bf{rt,rtp1});
        
        if rt==rtp1
            
            Grc=Grc+kron(speye(nzz),B{rt});
            
        end
        
        % Grc(abs(Grc)<1e-9)=0;
        
        [ir,jr,ss{iter}]=find(Grc);
        
        ii{iter}=offsetr+ir;
        
        jj{iter}=offsetc+jr;
        
        offsetc=c_range(end);
        
    end
    
    offsetr=r_range(end);
    
end

ii=cell2mat(ii);

jj=cell2mat(jj);

ss=cell2mat(ss);

G=sparse(ii,jj,ss,nt*nzz*h,nt*nzz*h,numel(ss));

clear ii jj ss A C B fv_lamba_bf

Tzz0=G\sparse(a);

retcode=1-all(isfinite(Tzz0));

Tzz=cell(1,h);

the_range=1:nt*nzz;

for rt=1:h
    
    Tzz{rt}=reshape(Tzz0(the_range),nt,nzz);
    
    if need_uncompress
        
        Tzz{rt}=Tzz{rt}*W;
        
    end
    
    if rt<h
        
        the_range=the_range+nt*nzz;
        
    end
    
end

end