function obj=assign_estimates(obj,params)
% assign_estimates -- assign parameters to the model typically under
% estimation
%
% Syntax
% -------
% ::
%
%   obj=assign_estimates(obj,params)
%
% Inputs
% -------
%
% - **obj** [rise_generic]: model object
%
% - **params** [vector]: values of parameters estimated or under estimation
% or under posterior simulation
%
% Outputs
% --------
%
% - **obj** [rise_generic]: model object
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

% this is the routine called during estimation for assigning new estimates
% to the model
if isempty(obj)
    obj=struct();
    return
end

if ~isempty(params)
    % this is general enough and includes measurement errors
    ids=obj.estimation_restrictions(:,1);
    inflate=obj.estimation_restrictions(:,2);
    obj.parameter_values(ids)=params(inflate);
end
