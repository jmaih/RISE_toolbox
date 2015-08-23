function x=unstransform_estimates(obj,x)
% unstransform_estimates -- sets the transformed parameters back normal
%
% Syntax
% -------
% ::
%
%   x=unstransform_estimates(obj,x)
%
% Inputs
% -------
%
% - **obj** [rise|dsge|rfvar|svar]: model object
%
% - **x** [n x 1 vector]: original parameters
%
% Outputs
% --------
%
% - **x** [n x 1 vector]: untransformed parameters
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

if isempty(obj)
    if nargout>1
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    end
    x=struct();
    return
end

% the transformations occur in setup_priors in function do_the_dirichlet
%------------------------------------------------------------------------
% Now we undo them
for id=1:numel(obj.estim_dirichlet)
    pos=obj.estim_dirichlet(id).location;
    ax_diag=obj.estim_dirichlet(id).ax_diag;
    x(pos)=utils.distrib.dirichlet_untransform(x(pos),ax_diag);
end

end