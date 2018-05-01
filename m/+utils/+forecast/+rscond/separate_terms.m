function [Ty,Te,Tsig,C]=separate_terms(T,sstate,state_cols,k,nshocks)
% SEPARATE_TERMS -- separate terms for conditional forecasting
%
% ::
%
%
%   [Ty,Te,Tsig,C]=SEPARATE_TERMS(T,sstate,ny,nx,nshocks,k,h)
%
% Args:
%
%    - **T** [cell array]: solution impact
%
%    - **sstate** [cell array]: steady state
%
%    - **state_cols** [vector]: location of state variables excluding shocks
%
%    - **k** [numeric]: number of forward steps (anticipation)
%
%    - **nshocks** [numeric]: number of shocks
%
% Returns:
%    :
%
%    - **Ty** [cell array]: includes square matrices for the impact of
%    endogenous variables
%
%    - **Te** [cell array]: includes matrices for the impact of shocks
%
%    - **Tsig** [cell array]: includes vectors summing the impact of
%    uncertainty and the trend
%
%    - **C** [cell array]: includes vectors summing the impact of
%    uncertainty, the trend and the steady state
%
% Note:
%
% Example:
%
%    See also:

h=numel(T);

ny=size(T{1},1);

nx=numel(state_cols);

C=cell(1,h);

Ty=cell(1,h);

Te=cell(1,h);

Tsig=cell(1,h);

Iy=eye(ny);

dim1Dist=ny;

dim2Dist=nshocks*ones(1,k+1);

for ireg=1:h
    
    Ti=T{ireg};
    
    Tx=Ti(:,1:nx);
    
    Ty{ireg}=add_zeros(Tx);
    
    % merge trend and the uncertainty impact
    %---------------------------------------
    Tsig{ireg}=Ti(:,nx+1);
    
    Tsig{ireg}=real(Tsig{ireg})+imag(Tsig{ireg});
    
    Te{ireg}=mat2cell(Ti(:,nx+2:end),dim1Dist,dim2Dist);
    
    C{ireg}=(Iy-Ty{ireg})*sstate{ireg}+Tsig{ireg};
    
end

    function Ty=add_zeros(Tx)
        
        Ty=zeros(ny);
        
        Ty(:,state_cols)=Tx;
        
    end

end