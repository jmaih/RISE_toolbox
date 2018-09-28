function [Tzz,retcode]=mytfqmr(A,B,C,fv_lamba_bf,W)

verbose=true;

need_uncompress = nargin==5;

h=numel(A);

[nt,nzz]=size(A{1});

N=nt*nzz*h;

rhs=zeros(nt,nzz,h);

smart_initialization=false;

x0=rhs;

for ii=1:h
    
    rhs(:,:,ii)=-A{ii};
    
    if smart_initialization
        
        x0(:,:,ii)=-B{ii}\A{ii};
        
    end
    
end

TOL=sqrt(eps);

MAXIT=min(N,500);

[Tzz0,retcode]=transpose_free_quasi_minimum_residual(...
    @left_hand_side,rhs(:),... % right hand side
    x0(:),... % initial guess
    TOL,... % tolerance level
    MAXIT,... % maximum number of iterations
    verbose);

max(abs(left_hand_side(Tzz0)-rhs(:)))

Tzz0=reshape(Tzz0,nt,nzz,h);

Tzz=cell(1,h);

for ii=1:h
    
    if need_uncompress
        
        Tzz{ii}=sparse(Tzz0(:,:,ii)*W);
        
    else
        
        Tzz{ii}=sparse(Tzz0(:,:,ii));
        
    end
    
end

    function r=left_hand_side(t)
        
        t=reshape(t,nt,nzz,h);
        
        r=zeros(nt,nzz,h);
        
        for rt=1:h
            
            r(:,:,rt)=B{rt}*t(:,:,rt);
            
            for rtp1=1:h
                
                r(:,:,rt)=r(:,:,rt)+fv_lamba_bf{rt,rtp1}*t(:,:,rtp1)*C{rt,rtp1};
                
            end
            
        end
        
        r=r(:);
        
    end

end
