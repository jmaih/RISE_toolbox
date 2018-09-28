function [sp,Z,Z_is_selector]=time_series_length(Z,T,H,Q,R)

sp=struct();

nH=size(H,3); sp.H=@(t)min(t,nH);

[m,~,nT]=size(T); sp.T=@(t)min(t,nT);

nQ=size(Q,3); sp.Q=@(t)min(t,nQ);

nR=size(R,3); sp.R=@(t)min(t,nR);

[~,ncZ,nZ]=size(Z);

Z_is_selector=ncZ~=m;

if Z_is_selector
    
    if islogical(Z)
        
        Z=find(Z);
        
    end
    
else
    
    sp.Z=@(t)min(t,nZ);
    
end

end
