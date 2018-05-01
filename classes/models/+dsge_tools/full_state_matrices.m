function [Aplus,A0,Aminus,T0]=full_state_matrices(siz,sm,T0)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

% builds state matrices for the system
% Aplus*X_{t+1}+A0*X_{t}+Aminus*X_{t-1}=0
% the variables are ordered in the order_var order i.e.
% [Stat.,Pred.,Both,Frwrd] 
if nargin<3
    
    T0=[];
    
end

h=siz.h;

nd=siz.nd;

ns=siz.ns;

np=siz.np;

nb=siz.nb;

nf=siz.nf;

Aplus=zeros(nd,nd,h,h);

A0=zeros(nd,nd,h);

Aminus=zeros(nd,nd,h);

for r0=1:h
    
    A0(:,:,r0)=sm.d0{r0};
    
    Aminus(:,ns+(1:np+nb),r0)=sm.dpb_minus{r0};
    
    for r1=1:h
        
        Aplus(:,ns+np+(1:nb+nf),r0,r1)=sm.dbf_plus{r0,r1};
        
    end
    
end

if nargout>3 && ~isempty(T0)
    
    tmp=T0;
    
    endo_nbr=size(tmp,1);
    
    T0=zeros(endo_nbr,endo_nbr,h);
    
    state_ids=ns+(1:np+nb);
    
    T0(:,state_ids,:)=tmp;
    
end

end