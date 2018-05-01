function varargout=u_higher_order_moments(siz)
% u_higher_order_moments -- computes moments Ekron(u,u,...,u)
%
% ::
%
%
%   [M1,M2,M3,...,M5]=u_higher_order_moments(siz)
%
% Args:
%
%    - siz : [struct] with fields
%      - np: number of predetermined
%      - nb: number of predetermined and forward-looking
%      - ne: number of shocks
%      - nz: (total) number of state variables (pred,bobth,sig,shocks)
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

shocks_id=siz.np+siz.nb+1+(1:siz.ne);

sig_id=siz.np+siz.nb+1;

o=nargout;

varargout=cell(1,o);

for iout=1:o

    if iout==2
    
        V=eye(siz.ne);
    
    else
        
        V=zeros([siz.ne*ones(1,iout),1]);
        
        if iout==4
        
            for ii=1:siz.ne
            
                V(ii,ii,ii,ii)=3;
                
                if ii<siz.ne
                
                    for jj=ii+1:siz.ne
                    
                        tmp=perms([ii,ii,jj,jj]);
                        
                        pos=sub2ind(siz.ne*ones(1,iout),tmp(:,1),tmp(:,2),tmp(:,3),tmp(:,4));
                        
                        V(pos)=1;
                    
                    end
                    
                end
                
            end
            
        end
        
    end
    
    nrows=siz.nz^iout;
    
    varargout{iout}=sparse(nrows,nrows);
    
    if rem(iout,2)==0
    
        vpos=collect_positions(shocks_id);
        
        sigcol=collect_positions(sig_id);
        
        varargout{iout}(vpos(:),sigcol(:))=V(:);
    
    end
    
end

    function pos=collect_positions(inpos)
    
        pos=false(siz.nz*ones(1,iout));
        
        locs=repmat({inpos},1,iout);
        
        pos(locs{:})=true;
    
    end

end

