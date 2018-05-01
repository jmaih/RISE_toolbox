function R=build_grid(R,vp,vectorized)
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

if nargin<3
    
    vectorized=true;
    
end

% see also mygrid

[rg,cg]=size(R);

rg_star=max(rg,1);

if vectorized
    
    Ip=(1:vp).';
    
    Ip=Ip(:,ones(rg_star,1));
    
    ii_=ones(vp,1)*(1:rg);
    
    R=[R(ii_(:),:),Ip(:)];

else
    
    Ip=repmat(transpose(1:vp),1,rg_star);
    
    G0=zeros(rg_star*vp,cg);
    
    for ii=1:rg
        
        G0((ii-1)*vp+1:ii*vp,:)=R(ii*ones(vp,1),:);
    
    end
    
    R=[G0,Ip(:)];
    
end

end