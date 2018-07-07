function estim_names=get_estimated_parameter_names(obj)
% INTERNAL FUNCTION
%

estim_names={obj.estimation.priors.name};

estim_names=update_estimated_parameter_names(obj,estim_names);

end

