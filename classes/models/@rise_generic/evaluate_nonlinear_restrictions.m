function g=evaluate_nonlinear_restrictions(obj)
% evaluate_nonlinear_restrictions - evaluates general restrictions 
%
% Syntax
% -------
% ::
%
%   g=evaluate_nonlinear_restrictions(obj)
%
% Inputs
% -------
%
% - **obj** [rise|dsge|rfvar|svar] : scalar of vector or RISE model objects
%
% Outputs
% --------
%
% - **g** [cell array] : value of restrictions for each object.
%
% More About
% ------------
%
% - The restrictions will be processed as g(x)<=0. But all the user has to
% do is to put zero where the restrictions are not violated!!!
%
% Examples
% ---------
%
% See also:

nobj=numel(obj);
if nobj==0
    g=struct();
else
    g=cell(1,nobj);
    for iobj=1:nobj
        if ~isempty(obj(iobj).options.estim_general_restrictions)
            if isempty(obj(iobj).general_restrictions_data)
                obj(iobj)=setup_general_restrictions(obj(iobj));
            end
            % take the absolute value: the violations are the non-zero
            % elements
            %-------------------------------------------------------------
            g{iobj}=abs(obj(iobj).general_restrictions_data(obj(iobj)));
        end
    end
end

end