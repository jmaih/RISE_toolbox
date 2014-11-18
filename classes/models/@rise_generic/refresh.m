function newobj=refresh(obj)
% refresh - refresh the options of an old object with a newer version of
%   the software
%
% Syntax
% -------
% ::
%
%   newobj=refresh(obj)
%
% Inputs
% -------
%
% - **obj** [rise|dsge|svar|rfvar]: model object
%
% Outputs
% --------
%
% - **obj** [rise|dsge|svar|rfvar]: refreshed model object
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

if isempty(obj)
    newobj=struct();
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