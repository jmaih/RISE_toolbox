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
% - **obj** [rise|dsge]: model object
%
% - **params** [vector]: values of parameters estimated or under estimation
% or under posterior simulation
%
% Outputs
% --------
%
% - **obj** [rise|dsge]: model object
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
    obj=assign_estimates@rise_generic(obj,params);
    % ensure the model will be re-solved no matter what
    %---------------------------------------------------
    obj.warrant_resolving = true;
end
