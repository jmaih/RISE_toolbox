function estim_names=get_estimated_parameter_names(obj)
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


estim_names={obj.estimation.priors.name};

estim_names=update_estimated_parameter_names(obj,estim_names);

end

