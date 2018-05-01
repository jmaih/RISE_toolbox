function obj=format_parameters(obj,ParameterizationArray)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:


param_names=obj.parameters.name;

if obj.is_dsge_var_model
    
    obj.dsge_prior_weight_id=find(strcmp('dsge_prior_weight',param_names),1);
    
end

if obj.is_hybrid_expectations_model
    % set all hybrid coefficients to 1
    pnames=obj.parameters.name;
    
    good=strncmp(pnames,'hbe_param_',10);
    
    pnames=pnames(good);
    
    np=numel(pnames);
    
    hbe_calibration=[pnames(:),num2cell(ones(np,1))];
    
    obj=setup_calibration(obj,hbe_calibration);
    
end

% baseline calibration
%---------------------
obj=setup_calibration(obj,ParameterizationArray.Calibration);

% set up priors
%--------------
obj=setup_priors(obj,ParameterizationArray.Priors,ParameterizationArray.error_control);

% measurement errors restrictions
%--------------------------------
 obj=setup_measurement_errors(obj);
 
% check that all the parameters in the model are in use
%------------------------------------------------------
not_in_use=param_names(~obj.parameters.is_in_use);

if ~isempty(not_in_use)
    
    warning([' the following parameters do not seem to affect model ',obj.filename,...
        '. You may want to discard them from model for tidiness'])
    
    disp(not_in_use)
    
end

end

