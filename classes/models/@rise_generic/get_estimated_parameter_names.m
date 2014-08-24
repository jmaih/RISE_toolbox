function estim_names=get_estimated_parameter_names(obj)

estim_names={obj.estimation.priors.name};
% estim_names=parser.param_name_to_valid_param_name({obj.estimation.priors.name});

estim_names=update_estimated_parameter_names(obj,estim_names);

end

