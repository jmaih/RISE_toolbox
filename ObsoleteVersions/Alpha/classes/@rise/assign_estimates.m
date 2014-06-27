function obj=assign_estimates(obj,params)
% this is the routine called during estimation for assigning new estimates
% to the model

if nargin>1 && ~isempty(params)
    % this is general enough and includes measurement errors
    ids=obj.estimation_restrictions(:,1);
    inflate=obj.estimation_restrictions(:,2);
    obj.parameter_values(ids)=params(inflate);
end

%{
if nargin>1 && ~isempty(params)
    % this is general enough and includes measurement errors
    ids=obj.estimation_restrictions(:,1);
    params=params(obj.estimation_restrictions(:,2));
    obj=obj.set_parameters({ids,params});
end
% obj=set(obj,'parameters',value)
%}
