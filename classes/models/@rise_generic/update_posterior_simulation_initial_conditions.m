function start=update_posterior_simulation_initial_conditions(obj,start,new_vcov,acceptance_rate)
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
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

if start.adapt_covariance
    npar=size(new_vcov,1);
    [R,p]=chol(new_vcov+obj.options.mcmc_diagcov_adjust_coef*eye(npar));
    if p==0
        start.Cs=transpose(R)/start.c;
    end
else
    start.c=utils.mcmc.retune_coefficient(start.c,obj.options.mcmc_target_range,acceptance_rate);
end
% go to next round
%-----------------
start.number_of_burns=0;
end
