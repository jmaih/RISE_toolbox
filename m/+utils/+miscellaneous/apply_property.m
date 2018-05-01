function newv=apply_property(type,refv,newv)
% apply_range -- makes to vectors have the same minimum and maximum
%
% ::
%
%
%   new=apply_range(type,refv,newv)
%
% Args:
%
%    - **refv** [vector]: reference vector
%
%    - **newv** [vector]: vector to modify
%
% Returns:
%    :
%
%    - **newv** [vector]: modified vector
%
% Note:
%
% Example:
%
%    See also:

switch lower(type)
    
    case 'range'
        
        if max(newv)==min(newv)
            
            ratio=.5;
            
        else
            
            ratio=(newv-min(newv))/(max(newv)-min(newv));
            
        end
        
        newv=min(refv)+ratio*(max(refv)-min(refv));
        
    case 'max'
        
        newv=newv-max(newv)+max(refv);
        
    case 'min'
        
        newv=newv-min(newv)+min(refv);
        
    otherwise
        
        error(['unknown property:: ',type])
        
end

end