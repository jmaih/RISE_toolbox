function obj=rehash(obj,force)
if nargin<2
    force=false;
end
% push the parameters in the images into their respective objects
% rebuild the parameter object
if ~obj(1).estimation_under_way||force
    nobj=numel(obj);
    number_of_regimes=obj(1).NumberOfRegimes;
    for ii=1:nobj
        obj(ii).parameters=cell2object(obj(ii).parameters_image,'rise_param');
        if ~obj(ii).is_optimal_policy_model && ~obj(ii).is_sticky_information_model
            % I wanted to push the steady state right way but sometimes the steady
            % state is not known before the model is totally solved as in the case
            % of optimal policy... and in particular, when the number of endogenous
            % variables changes after the solution
            endo_nbr=obj(ii).NumberOfEndogenous(1);
            tmp=mat2cell(obj(ii).steady_state_and_balanced_growth_path,ones(2*endo_nbr,1),number_of_regimes);
            tmp0=tmp(1:endo_nbr,:);
            % this thing below works well... even after checking many
            % times :-)
            [obj(ii).varendo(:).det_steady_state]=(tmp0{:});
            tmp0=tmp(endo_nbr+1:end,:);
            [obj(ii).varendo(:).balanced_growth]=(tmp0{:});
        end        
    end
end

end