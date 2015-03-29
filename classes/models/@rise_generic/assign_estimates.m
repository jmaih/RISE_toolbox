function obj=assign_estimates(obj,params)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
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
    % ensure the model will be re-solved no matter what
    %---------------------------------------------------
    if ~obj.estimation_under_way
        obj=set(obj,'parameters',{});
    end
end
