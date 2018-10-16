function obj=setup_endogenous_priors(obj,fh)
% INTERNAL FUNCTION
%

[obj.estim_endogenous_priors_data,...
obj.estimation.endogenous_priors,...
obj.options.estim_endogenous_priors]=utils.prior.setup_endogenous_priors_engine(obj.options.prior_trunc,fh);

end
