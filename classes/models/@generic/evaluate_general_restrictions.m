function g=evaluate_general_restrictions(obj)
% EVALUATE_GENERAL_RESTRICTIONS - evaluates general restrictions
%
% ::
%
%
%   g=EVALUATE_GENERAL_RESTRICTIONS(obj)
%
% Args:
%
%    - **obj** [rise|dsge|rfvar|svar] : scalar of vector or RISE model objects
%
% Returns:
%    :
%
%    - **g** [cell array] : value of restrictions for each object.
%
% Note:
%
%    - The restrictions will be processed as g(x)<=0. But all the user has to
%    do is to put zero where the restrictions are not violated!!!
%
% Example:
%
%    See also:

nobj=numel(obj);

if nobj==0
    
    g=cell(0,4);
    
else
    
    g=cell(1,nobj);
    
    for iobj=1:nobj
        
        if ~isempty(obj(iobj).options.estim_general_restrictions)
            % take the absolute value: the violations are the non-zero
            % elements
            %-------------------------------------------------------------
            g{iobj}=abs(obj(iobj).general_restrictions_data(obj(iobj)));
            
        end
        
    end
    
end

end