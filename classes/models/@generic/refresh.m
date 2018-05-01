function newobj=refresh(obj)
% REFRESH - refresh the options of an old object with a newer version of
%   the software
%
% ::
%
%
%   newobj=REFRESH(obj)
%
% Args:
%
%    - **obj** [rise|dsge|svar|rfvar]: model object
%
% Returns:
%    :
%
%    - **newobj** [rise|dsge|svar|rfvar]: refreshed model object
%
% Note:
%
% Example:
%
%    See also:

if isempty(obj)
    
    newobj=cell(0,4);
    
else
    
    newobj=obj;
    
    newobj.options=[];
    
    newobj=set(newobj,'initialize');
    
    default_options=newobj.options;
    
    old_options=obj.options;
    
    newfields=fieldnames(default_options);
    
    oldfields=fieldnames(old_options);
    
    % take the old fields into the new ones but on the new terms
    %-----------------------------------------------------------
    settable=ismember(newfields,oldfields);
    
    for ifield=1:numel(newfields)
        
        if settable(ifield)
            
            ff=newfields{ifield};
            
            default_options.(ff)=old_options.(ff);
            
        end
        
    end
    
    newobj.options=default_options;
    
end

end