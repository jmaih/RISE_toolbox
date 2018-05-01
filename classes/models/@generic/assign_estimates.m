function obj=assign_estimates(obj,params)
% assign_estimates -- assign parameters to the model typically under
% estimation
%
% ::
%
%
%   obj=assign_estimates(obj,params)
%
% Args:
%
%    - **obj** [rise|dsge|svar|rfvar]: model object
%
%    - **params** [vector]: values of parameters estimated or under estimation
%    or under posterior simulation
%
% Returns:
%    :
%
%    - **obj** [generic]: model object
%
% Note:
%
% Example:
%
%    See also:

% this is the routine called during estimation for assigning new estimates
% to the model
if isempty(obj)
    obj=struct();
    return
end

if ~isempty(params)
    ids=obj.estimation_restrictions(:,1);
    inflate=obj.estimation_restrictions(:,2);
    % this is general enough and includes measurement errors
    M=obj.parameter_values;
    % push estimates
    %---------------
    M(ids)=params(inflate);
    % add the derived parameters
    %----------------------------
    if ~isempty(obj.routines.derived_parameters)
        linrest=obj.routines.derived_parameters;
        for irest=1:size(linrest,1)
            row=linrest{irest,1};
            cols=linrest{irest,2};
            M(row,cols)=linrest{irest,3}(M);
        end
    end
    % write back to the object
    %--------------------------
    obj.parameter_values=M;
end
