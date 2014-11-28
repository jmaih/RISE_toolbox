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
% - **g** [cell array] : value of restrictions for each object
%
% More About
% ------------
%
% - The restrictions are expected to be of the form g(x)<=0. This must be
% ensured by the user who writes the restriction function to be passed to
% RISE.
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
            g{iobj}=obj(iobj).general_restrictions_data(obj(iobj));
        end
    end
end

end