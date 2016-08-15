function [pos,siz,shock_horizon]=rehash_topology(obj,structural_matrices)

pos=obj.locations.before_solve;

siz=obj.siz.before_solve;

% collect the sizes
%------------------

siz.h=size(structural_matrices.dv,2);

siz.nT=siz.ns+siz.np+siz.nb+siz.nf;

shock_horizon=max(obj.exogenous.shock_horizon(:));

siz.nz=siz.np+siz.nb+1+siz.ne*(1+shock_horizon);

siz.nd=size(structural_matrices.dv{1,1},1); % number of equations

if siz.ne
    
    pos.z.e_plus=pos.z.e_0(end)+(1:shock_horizon*siz.ne);
    
else
    
    pos.z.e_plus=(1:shock_horizon*siz.ne);
    
end


end