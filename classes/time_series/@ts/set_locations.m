function locs=set_locations(dnx,dn)
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

% find the locations of dnx in dn. It is assumed that dn has
% consecutive dates.
if ~is_consecutive(dn)
    
    error('the numbers in the second argument should be consecutive')
    
end

if is_consecutive(dnx)
    
    first=find(dnx(1)==dn);
    
    last=find(dnx(end)==dn);
    
    if isempty(first)||isempty(last)
        
        error('all date numbers in the first should be part of the second')
        
    end
    
    locs=first:last;% <--locs=first+(0:numel(dnx)-1);
    
else
    
    locs=nan(1,numel(dnx));
    
    for irow=1:numel(dnx)
        
        tmp=find(dnx(irow)==dn);
        
        if isempty(tmp)
            
            error('all date numbers in the first should be part of the second')
            
        end
        
        locs(irow)=tmp;
        
    end
    
end

    function flag=is_consecutive(d)
        
        flag=true;
        
        if numel(d)>1
            
            test=d(2:end)-d(1:end-1);
            
            flag=all(test(:)-1<sqrt(eps));
            
        end
        
    end

end
