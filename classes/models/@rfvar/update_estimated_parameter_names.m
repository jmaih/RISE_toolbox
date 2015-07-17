function estim_names=update_estimated_parameter_names(obj,estim_names)
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


is_constant_parameter_var=obj.markov_chains.regimes_number==1 && ...
        obj.options.vp_analytical_post_mode;
% if constant paramater and analytical solution, keep only the ai and ci
% parameters
if is_constant_parameter_var
    [lag_names]=vartools.select_parameter_type(estim_names,'lag_coef');
    [det_names]=vartools.select_parameter_type(estim_names,'det_coef');
    estim_names=[lag_names,det_names];
end

end

