function x=unstransform_parameters(obj,x)
% UNTRANSFORM_PARAMETERS -- sets the transformed parameters back normal
%
% ::
%
%
%   x=UNTRANSFORM_PARAMETERS(obj,x)
%
% Args:
%
%    - **obj** [rise|dsge|rfvar|svar]: model object
%
%    - **x** [n x 1 vector]: original parameters
%
% Returns:
%    :
%
%    - **x** [n x 1 vector]: untransformed parameters
%
% Note:
%
% Example:
%
%    See also: TRANSFORM_PARAMETERS

if isempty(obj)
    
    if nargout>1
    
        error([mfilename,':: when the object is emtpy, nargout must be at most 1'])
    
    end
    
    x=struct();
    
    return

end

linear_restricts=obj(1).linear_restrictions_data;
% expand x before doing anything
%-------------------------------
x=linear_restricts.a_func(x);

% the transformations occur in setup_priors in function do_the_dirichlet
%------------------------------------------------------------------------
% Now we undo them
for id=1:numel(obj.estim_priors_data.estim_dirichlet)
    
    pos=obj.estim_priors_data.estim_dirichlet(id).location;
    
    sum_aij=obj.estim_priors_data.estim_dirichlet(id).sum_aij;
    
    x(pos)=utils.distrib.dirichlet_untransform(x(pos),sum_aij);
    
end

end