function Tz_sig=solve_perturbation_impact(Tz_sig,A0sig,dbf_plus,dt_t,siz,pos)

if any(dt_t(:))% then automatically h>1
    % use a qr decomposition to solve a small system. Given the structure
    % of the system, it is enough to precondition it.
    %-----------------------------------------------
    for r0=1:siz.h
        
        A0_sig_i=A0sig(:,:,r0)\eye(siz.nd);
        
        for r1=1:siz.h
            
            dbf_plus{r0,r1}=A0_sig_i*dbf_plus{r0,r1};
            
        end
        
        dt_t(:,r0)=A0_sig_i*dt_t(:,r0);
        
    end
    % now we solve the system sum(A+*Tz_sig(+)+Tz_sig+dt_t=0 first for
    % variables p,b,f and then for variables s
    clear A0sig
    
    % solve the small system without static variables
    %------------------------------------------------
    % the direct solution implemented below is not efficient in very large
    % systems...
    npbf=siz.np+siz.nb+siz.nf;
    
    A=zeros(npbf*siz.h);
    
    for r0=1:siz.h
        
        row_=(r0-1)*npbf+1:r0*npbf;
        
        for r1=1:siz.h
            
            col_=(r1-1)*npbf+1:r1*npbf;
            
            A(row_,col_)=[zeros(npbf,siz.np),dbf_plus{r0,r1}(siz.ns+1:end,:)];
            
            if r0==r1
                
                A(row_,col_)=A(row_,col_)+eye(npbf);
                
            end
            
        end
        
    end
    
    Tz_sig_PBF=uminus(dt_t(siz.ns+1:end,:));
    
    Tz_sig_PBF=A\Tz_sig_PBF(:);
    
    Tz_sig(siz.ns+1:end,:)=reshape(Tz_sig_PBF,npbf,siz.h);
    
    % solve the static variables
    %---------------------------
    for r0=1:siz.h
        
        Tz_sig(1:siz.ns,r0)=dt_t(1:siz.ns,r0);
        
        for r1=1:siz.h
            
            Tz_sig(1:siz.ns,r0)=Tz_sig(1:siz.ns,r0)+...
                dbf_plus{r0,r1}(1:siz.ns,:)*Tz_sig(pos.t.bf,r1);
            
        end
        
        Tz_sig(1:siz.ns,r0)=uminus(Tz_sig(1:siz.ns,r0));
        
    end
    
end

end
