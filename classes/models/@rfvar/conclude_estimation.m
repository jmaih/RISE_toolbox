function [log_post,log_lik,log_prior,Incr,retcode,obj,x11]=conclude_estimation(obj,x1)
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

if obj.markov_chains.regimes_number==1 && obj.options.estim_analytical_post_mode
    % then the vector is not full
    %----------------------------
    vdata=obj.constant_var_data.vdata;
    
    vcov=obj.constant_var_data.vcov;
    
    x11=vartools.build_parameter_vector(vdata,x1,vcov);
    
else
    
    x11=x1;
    
end

[log_post,log_lik,log_prior,Incr,retcode,obj,x11]=conclude_estimation@svar(obj,x11);

end